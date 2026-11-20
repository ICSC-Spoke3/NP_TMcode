/* Copyright (C) 2025   INAF - Osservatorio Astronomico di Cagliari

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   A copy of the GNU General Public License is distributed along with
   this program in the COPYING file. If not, see: <https://www.gnu.org/licenses/>.
 */

/**
 * \file cublas_calls.cpp
 *
 * \brief Implementation of the interface with CUBLAS libraries.
 *
 * The matrix is inverted by LU decomposition with cuSOLVER (`zgetrf`, then
 * `zgetrs` on the identity); cuBLAS runs the refinement.
 *
 * Each thread keeps its GPU resources (stream, library handles and device
 * buffers) between inversions, because creating them for every matrix costs
 * more than the inversion of a small matrix. They are released when a GPU
 * call fails, and left to the operating system at the end of the program.
 * \ingroup linalg
 */

#ifdef USE_CUBLAS

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <vector>
// In clang's OpenMP offload compilation for NVIDIA GPUs, the CUDA headers
// define __noinline__ as a macro, which breaks the libstdc++ headers (GCC 12
// and later) that use __attribute__((__noinline__)). Include those first,
// and remove the macro after the CUDA headers, unless compiling CUDA code.
#include <memory>
#include <string>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#if defined(__noinline__) && !defined(__CUDACC__)
#undef __noinline__
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_CUBLAS_CALLS_H_
#include "../include/cublas_calls.h"
#endif

#ifdef USE_ILP64
//! \brief cuBLAS complex matrix-matrix product (64-bit integer interface).
#define CUZGEMM cublasZgemm_64
//! \brief cuBLAS complex vector update y = a * x + y (64-bit integer interface).
#define CUZAXPY cublasZaxpy_64
//! \brief cuBLAS index of the complex element with maximum modulus (64-bit integer interface).
#define CUIZAMAX cublasIzamax_64
//! \brief Integer type of the cuBLAS interface.
#define CUINT int64_t
#else
//! \brief cuBLAS complex matrix-matrix product.
#define CUZGEMM cublasZgemm
//! \brief cuBLAS complex vector update y = a * x + y.
#define CUZAXPY cublasZaxpy
//! \brief cuBLAS index of the complex element with maximum modulus.
#define CUIZAMAX cublasIzamax
//! \brief Integer type of the cuBLAS interface.
#define CUINT int
#endif

//! \cond
namespace {
  std::string cuda_failure(const char *what, int line) {
    return std::string(what) + " (cublas_calls.cpp, line " + std::to_string(line) + ")";
  }
}
//! \endcond

//! \brief Throw a DeviceError if a CUDA runtime call fails.
#define cudacall(call)							\
  do {									\
    cudaError_t err = (call);						\
    if (cudaSuccess != err) {						\
      throw DeviceError(cuda_failure(cudaGetErrorString(err), __LINE__)); \
    }									\
  } while (0)

//! \brief Throw a DeviceError if a cuBLAS call fails.
#define cublascall(call)						\
  do {									\
    cublasStatus_t status = (call);					\
    if (CUBLAS_STATUS_SUCCESS != status) {				\
      throw DeviceError(cuda_failure(					\
	("cuBLAS status " + std::to_string((int)status)).c_str(), __LINE__)); \
    }									\
  } while (0)

//! \brief Throw a DeviceError if a cuSOLVER call fails.
#define cusolvercall(call)						\
  do {									\
    cusolverStatus_t status = (call);					\
    if (CUSOLVER_STATUS_SUCCESS != status) {				\
      throw DeviceError(cuda_failure(					\
	("cuSOLVER status " + std::to_string((int)status)).c_str(), __LINE__)); \
    }									\
  } while (0)

using namespace std;

//! \cond
namespace {
  // Device memory that grows on demand and is kept between inversions.
  class DeviceBuffer {
  public:
    void *ptr;
    size_t bytes;
    DeviceBuffer() : ptr(NULL), bytes(0) { }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    void release() {
      if (ptr != NULL) cudaFree(ptr);
      ptr = NULL;
      bytes = 0;
    }
    // Make room for size bytes; false if the GPU has no room.
    bool reserve(size_t size) {
      if (size <= bytes) return true;
      release();
      if (cudaMalloc(&ptr, size) != cudaSuccess) {
	cudaGetLastError(); // clear the allocation error
	ptr = NULL;
	return false;
      }
      bytes = size;
      return true;
    }
    // As reserve(), with a DeviceError if the GPU has no room.
    void require(size_t size, const char *what) {
      if (!reserve(size)) {
	throw DeviceError(string("not enough GPU memory for ") + what);
      }
    }
    template <class T> T* as() { return (T *)ptr; }
  };

  // GPU resources of one thread on one device.
  class GpuContext {
  public:
    int device;
    bool ready;
    cudaStream_t stream;
    cublasHandle_t blas;
    cusolverDnHandle_t solver;
    cusolverDnParams_t params;
    std::vector<char> host_work; // host workspace of cuSOLVER
    DeviceBuffer matrix, inverse, residual, refine, one, work, pivots, info;
    GpuContext() : device(-1), ready(false), stream(NULL), blas(NULL), solver(NULL), params(NULL) { }
    GpuContext(const GpuContext&) = delete;
    GpuContext& operator=(const GpuContext&) = delete;

    // Create the stream and the handles on device_id, unless they exist.
    void setup(int device_id) {
      if (ready && device == device_id) return;
      reset();
      cudacall(cudaSetDevice(device_id));
      device = device_id;
      cudacall(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
      cublascall(cublasCreate_v2(&blas));
      cublascall(cublasSetStream_v2(blas, stream));
      cusolvercall(cusolverDnCreate(&solver));
      cusolvercall(cusolverDnSetStream(solver, stream));
      cusolvercall(cusolverDnCreateParams(&params));
      ready = true;
    }

    // Release everything (after a failure, the state of the GPU is unknown).
    void reset() {
      if (device >= 0) cudaSetDevice(device);
      matrix.release(); inverse.release(); residual.release(); refine.release();
      one.release(); work.release(); pivots.release(); info.release();
      if (params != NULL) cusolverDnDestroyParams(params);
      if (solver != NULL) cusolverDnDestroy(solver);
      params = NULL;
      solver = NULL;
      host_work.clear();
      host_work.shrink_to_fit();
      if (blas != NULL) cublasDestroy_v2(blas);
      if (stream != NULL) cudaStreamDestroy(stream);
      blas = NULL;
      stream = NULL;
      ready = false;
      device = -1;
      cudaGetLastError(); // clear any error left by the failure
    }

    void sync() { cudacall(cudaStreamSynchronize(stream)); }
  };

  // One context per thread. It is never destroyed: the CUDA runtime may
  // already be unloaded when thread-local objects are destroyed at exit.
  GpuContext& thread_context() {
    thread_local GpuContext *context = new GpuContext();
    return *context;
  }

  // Copy the inverse to the host matrix. A failure during the copy may
  // leave the host matrix partly overwritten, so it is not a DeviceError
  // (which would let the caller invert the matrix again on the CPU).
  void copy_to_host(GpuContext &ctx, cuDoubleComplex *a, size_t bytes) {
    cudaError_t result = cudaMemcpyAsync(a, ctx.inverse.ptr, bytes, cudaMemcpyDeviceToHost, ctx.stream);
    if (result == cudaSuccess) result = cudaStreamSynchronize(ctx.stream);
    if (result != cudaSuccess) {
      string what = cudaGetErrorString(result);
      ctx.reset();
      throw NptmError("the GPU failed while copying the inverse matrix to the host (" + what + ")", NPTM_EXIT_DEVICE);
    }
  }

  void log_message(const RuntimeSettings& rs, const string& message, int level = LOG_INFO) {
    if (rs.logger != NULL) rs.logger->log(message, level);
  }

  double modulus(const cuDoubleComplex &value) {
    return std::hypot(cuCreal(value), cuCimag(value));
  }

  // Invert the m x m matrix in ctx.matrix into ctx.inverse (both
  // column-major for cuBLAS; the inverse of the transpose is the transpose
  // of the inverse, so row-major data give the right result). The content
  // of ctx.matrix is destroyed. Returns the info code of the decomposition
  // (> 0 for a singular matrix).
  int invert_on_device(GpuContext &ctx, CUINT m, const cuDoubleComplex &cu_one) {
    int h_info = 0;
    ctx.info.require(sizeof(int), "the inversion status");
    const size_t matrix_bytes = (size_t)m * (size_t)m * sizeof(cuDoubleComplex);
    // LU decomposition, then solution of A X = I.
    size_t device_work = 0, host_work = 0;
    cusolvercall(cusolverDnXgetrf_bufferSize(
      ctx.solver, ctx.params, (int64_t)m, (int64_t)m, CUDA_C_64F, ctx.matrix.ptr, (int64_t)m,
      CUDA_C_64F, &device_work, &host_work
    ));
    ctx.work.require(device_work > 0 ? device_work : 1, "the LU workspace");
    if (ctx.host_work.size() < host_work) ctx.host_work.resize(host_work);
    ctx.pivots.require((size_t)m * sizeof(int64_t), "the pivot indices");
    cusolvercall(cusolverDnXgetrf(
      ctx.solver, ctx.params, (int64_t)m, (int64_t)m, CUDA_C_64F, ctx.matrix.ptr, (int64_t)m,
      ctx.pivots.as<int64_t>(), CUDA_C_64F, ctx.work.ptr, device_work,
      host_work > 0 ? ctx.host_work.data() : NULL, host_work, ctx.info.as<int>()
    ));
    cudacall(cudaMemcpyAsync(&h_info, ctx.info.ptr, sizeof(int), cudaMemcpyDeviceToHost, ctx.stream));
    ctx.sync();
    if (h_info < 0) {
      throw NptmError("cusolverDnXgetrf was called with illegal argument " + to_string(-h_info));
    }
    if (h_info > 0) return h_info;
    // The identity: zeros, then ones on the diagonal (stride m + 1).
    cudacall(cudaMemsetAsync(ctx.inverse.ptr, 0, matrix_bytes, ctx.stream));
    cublascall(CUZAXPY(ctx.blas, m, &cu_one, ctx.one.as<cuDoubleComplex>(), 0, ctx.inverse.as<cuDoubleComplex>(), m + 1));
    cusolvercall(cusolverDnXgetrs(
      ctx.solver, ctx.params, CUBLAS_OP_N, (int64_t)m, (int64_t)m, CUDA_C_64F, ctx.matrix.ptr,
      (int64_t)m, ctx.pivots.as<int64_t>(), CUDA_C_64F, ctx.inverse.ptr, (int64_t)m, ctx.info.as<int>()
    ));
    cudacall(cudaMemcpyAsync(&h_info, ctx.info.ptr, sizeof(int), cudaMemcpyDeviceToHost, ctx.stream));
    ctx.sync();
    if (h_info < 0) {
      throw NptmError("cusolverDnXgetrs was called with illegal argument " + to_string(-h_info));
    }
    return 0;
  }

  // Newton-Schulz refinement of ctx.inverse, with the original matrix in
  // ctx.matrix. Returns 0 if the refined inverse is to be kept, 1 if the
  // unrefined one (already on the host) is to be kept.
  int refine_on_device(GpuContext &ctx, CUINT m, const RuntimeSettings& rs) {
    const CUINT mm = m * m;
    char buffer[128];
    int ref_err = 0;
    cuDoubleComplex cu_mone = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex cu_one = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex cu_zero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex *d_a = ctx.matrix.as<cuDoubleComplex>();
    cuDoubleComplex *d_c = ctx.inverse.as<cuDoubleComplex>();
    cuDoubleComplex *d_residual = ctx.residual.as<cuDoubleComplex>();
    cuDoubleComplex *d_refine = ctx.refine.as<cuDoubleComplex>();
    // Detect the maximum value of the inverse matrix.
    double oldmax = 2.0e+16, curmax = 1.0e+16;
    CUINT maxindex = 1;
    cuDoubleComplex cublasmax;
    cublascall(CUIZAMAX(ctx.blas, mm, d_c, 1, &maxindex));
    cudacall(cudaMemcpyAsync(&cublasmax, d_c + maxindex - 1, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost, ctx.stream));
    ctx.sync();
    curmax = modulus(cublasmax);
    snprintf(buffer, sizeof(buffer), "largest matrix value has modulus %.5le.\n", curmax);
    log_message(rs, buffer);
    // Iterative refinement loop
    int max_ref_iters = rs.max_ref_iters;
    for (int ri = 0; ri < max_ref_iters; ri++) {
      oldmax = curmax;
      // multiply minus the original matrix times the inverse matrix
      // NOTE: factors in zgemm are swapped because zgemm is designed for column-major
      // Fortran-style arrays, whereas our arrays are C-style row-major.
      cublascall(CUZGEMM(ctx.blas, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_mone, d_c, m, d_a, m, &cu_zero, d_residual, m));
      // add the identity to the product
      cublascall(CUZAXPY(ctx.blas, m, &cu_one, ctx.one.as<cuDoubleComplex>(), 0, d_residual, m + 1));
      // find the maximum absolute value of the residual
      cublascall(CUIZAMAX(ctx.blas, mm, d_residual, 1, &maxindex));
      // transfer the maximum value to the host
      cudacall(cudaMemcpyAsync(&cublasmax, d_residual + maxindex - 1, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost, ctx.stream));
      ctx.sync();
      curmax = modulus(cublasmax);
      snprintf(buffer, sizeof(buffer), "iteration %d has residue %.5le; target residue is %.5le.\n", ri, curmax, rs.accuracy_goal);
      log_message(rs, buffer, LOG_DEBG);
      if (curmax < 0.5) { // Safe conditions for Newton-Schultz iteration.
	if (curmax < 0.95 * oldmax) { // Newton-Schultz iteration is improving and can proceed.
	  // multiply the inverse times the residual, add to the initial inverse
	  cublascall(CUZGEMM(ctx.blas, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_one, d_residual, m, d_c, m, &cu_zero, d_refine, m));
	  // add to the initial inverse
	  cublascall(CUZAXPY(ctx.blas, mm, &cu_one, d_refine, 1, d_c, 1));
	  if (curmax < rs.accuracy_goal) {
	    log_message(rs, "good news - optimal convergence achieved. Stopping.\n", LOG_DEBG);
	    ref_err = 0;
	    break; // ri for
	  }
	} else {
	  if (curmax > 0.1) {
	    snprintf(buffer, sizeof(buffer), "iteration %d achieved limiting residue %.5le. Cannot reach goal. Reverting.\n", ri, curmax);
	    log_message(rs, buffer);
	    ref_err = 1;
	  } else {
	    snprintf(buffer, sizeof(buffer), "iteration %d achieved limiting residue %.5le. Stopping.\n", ri, curmax);
	    log_message(rs, buffer, LOG_WARN);
	  }
	  break; // ri for
	}
      } else { // curmax > 0.5. Newton-Schultz iteration is dangerous.
	if (curmax < oldmax) {
	  snprintf(buffer, sizeof(buffer), "iteration %d has residue %.5le. Iterating is dangerous. Stopping.\n", ri, curmax);
	  log_message(rs, buffer, LOG_WARN);
	} else {
	  ref_err = 1;
	  snprintf(buffer, sizeof(buffer), "iteration %d has residue %.5le. Reverting to unrefined and stopping.\n", ri, curmax);
	  log_message(rs, buffer);
	}
	break; // ri for
      }
    } // end of ri for
    return ref_err;
  }
}
//! \endcond

void cublas_zinvert(dcomplex **mat, np_int n, int &jer, int device_id, const RuntimeSettings& rs) {
  // The host matrix changes only when the inversion has succeeded, or holds
  // the unrefined inverse during the refinement, so that a caller can invert
  // it in another way if a DeviceError is thrown before.
  jer = 0;
  if (rs.invert_mode != RuntimeSettings::INV_MODE_LU) {
    throw NptmError("the cuBLAS solver only implements LU inversion"); // check_inversion_mode() rejects the others
  }
#ifndef USE_ILP64
  if (n > (np_int)numeric_limits<int>::max()) {
    throw DeviceError("the matrix is too large for the 32-bit cuBLAS interface");
  }
#endif
  const CUINT m = (CUINT)n;
  const size_t matrix_bytes = (size_t)m * (size_t)m * sizeof(cuDoubleComplex);
  cuDoubleComplex *a = (cuDoubleComplex *)&(mat[0][0]); // pointer to first element on host
  const cuDoubleComplex cu_one = make_cuDoubleComplex(1.0, 0.0);
  GpuContext &ctx = thread_context();
  try {
    ctx.setup(device_id);
    cudacall(cudaSetDevice(device_id));
    ctx.matrix.require(matrix_bytes, "the matrix");
    ctx.inverse.require(matrix_bytes, "the inverse matrix");
    if (ctx.one.bytes < sizeof(cuDoubleComplex)) {
      ctx.one.require(sizeof(cuDoubleComplex), "a constant");
      cudacall(cudaMemcpyAsync(ctx.one.ptr, &cu_one, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice, ctx.stream));
    }
    cudacall(cudaMemcpyAsync(ctx.matrix.ptr, a, matrix_bytes, cudaMemcpyHostToDevice, ctx.stream)); // copy a -> d_a
    int info = invert_on_device(ctx, m, cu_one);
    if (info > 0) {
      log_message(rs, "the cuBLAS inversion returned code " + to_string(info) + " (singular matrix).\n", LOG_WARN);
      jer = info;
      return; // the host matrix is unchanged
    }
  } catch (...) {
    ctx.reset();
    throw;
  }
  int ref_err = 0;
  bool host_has_inverse = false;
  if (rs.use_refinement) {
    // The refinement buffers are allocated before the host matrix changes:
    // if the GPU has no room, the inverse is not refined.
    if (!ctx.residual.reserve(matrix_bytes) || !ctx.refine.reserve(matrix_bytes)) {
      ctx.residual.release();
      ctx.refine.release();
      log_message(rs, "not enough GPU memory for the refinement: the inverse is not refined.\n", LOG_WARN);
    } else {
      try {
	// The original matrix again in d_a, and the unrefined inverse on the
	// host as a back-up if the refinement fails.
	cudacall(cudaMemcpyAsync(ctx.matrix.ptr, a, matrix_bytes, cudaMemcpyHostToDevice, ctx.stream));
	ctx.sync();
      } catch (...) {
	ctx.reset();
	throw; // the host matrix is still unchanged
      }
      copy_to_host(ctx, a, matrix_bytes);
      host_has_inverse = true;
      // From here on, the host holds the unrefined inverse: a GPU failure stops the refinement.
      try {
	ref_err = refine_on_device(ctx, m, rs);
      } catch (const DeviceError &ex) {
	log_message(
	  rs, string("the GPU failed during the refinement (") + ex.what() + "): the inverse is not refined.\n",
	  LOG_WARN
	);
	ctx.reset();
	return; // keep the unrefined inverse on the host
      }
    }
  }
  if (ref_err == 0) { // Copy the (refined) inverse back to the host.
    try {
      ctx.sync(); // a failure here leaves the host matrix as it was
    } catch (...) {
      ctx.reset();
      if (!host_has_inverse) throw;
      log_message(rs, "the GPU failed after the refinement: the inverse is not refined.\n", LOG_WARN);
      return; // keep the unrefined inverse on the host
    }
    copy_to_host(ctx, a, matrix_bytes);
  }
}


#endif
