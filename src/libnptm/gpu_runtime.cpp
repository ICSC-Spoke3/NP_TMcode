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

/*! \file gpu_runtime.cpp
 *
 * \brief Implementation of the GPU runtime calls, for CUDA or HIP.
 * \ingroup support
 */

#include <cstddef>
#include <memory>
#include <string>

// The runtime of the build: the one MAGMA was built on, or CUDA for cuBLAS.
// In clang's OpenMP offload compilation for NVIDIA GPUs, the CUDA headers
// define __noinline__ as a macro, which breaks the libstdc++ headers (GCC 12
// and later): the C++ headers are included first, and the macro is removed.
#if defined(USE_MAGMA)
#include <magma_v2.h>
#if defined(MAGMA_HAVE_HIP)
#define NPTM_GPU_HIP
#else
#define NPTM_GPU_CUDA
#endif
#elif defined(USE_CUBLAS)
#define NPTM_GPU_CUDA
#endif

#if defined(NPTM_GPU_HIP)
#include <hip/hip_runtime_api.h>
#elif defined(NPTM_GPU_CUDA)
#include <cuda_runtime.h>
#endif

#if defined(__noinline__) && !defined(__CUDACC__)
#undef __noinline__
#endif

#ifndef INCLUDE_GPU_RUNTIME_H_
#include "../include/gpu_runtime.h"
#endif

const char* gpu_runtime_name() {
#if defined(NPTM_GPU_HIP)
  return "HIP";
#elif defined(NPTM_GPU_CUDA)
  return "CUDA";
#else
  return "none";
#endif
}

int gpu_device_count(std::string& error) {
  int count = 0;
#if defined(NPTM_GPU_HIP)
  hipError_t result = hipGetDeviceCount(&count);
  if (result != hipSuccess) {
    error = std::string("HIP: ") + hipGetErrorString(result);
    return -1;
  }
#elif defined(NPTM_GPU_CUDA)
  cudaError_t result = cudaGetDeviceCount(&count);
  if (result != cudaSuccess) {
    error = std::string("CUDA: ") + cudaGetErrorString(result);
    return -1;
  }
#else
  error = "the build uses no GPU runtime";
  return -1;
#endif
  return count;
}

long gpu_free_memory(int device) {
  long result = -1;
#if defined(NPTM_GPU_HIP) || defined(NPTM_GPU_CUDA)
  size_t free_bytes = 0, total_bytes = 0;
  int current = 0;
#endif
#if defined(NPTM_GPU_HIP)
  if (hipGetDevice(&current) == hipSuccess && hipSetDevice(device) == hipSuccess) {
    if (hipMemGetInfo(&free_bytes, &total_bytes) == hipSuccess) result = (long)free_bytes;
    hipSetDevice(current);
  }
#elif defined(NPTM_GPU_CUDA)
  if (cudaGetDevice(&current) == cudaSuccess && cudaSetDevice(device) == cudaSuccess) {
    if (cudaMemGetInfo(&free_bytes, &total_bytes) == cudaSuccess) result = (long)free_bytes;
    cudaSetDevice(current);
  }
#else
  (void)device;
#endif
  return result;
}

bool gpu_set_stack_limit(size_t bytes) {
#if defined(NPTM_GPU_HIP)
  return hipDeviceSetLimit(hipLimitStackSize, bytes) == hipSuccess;
#elif defined(NPTM_GPU_CUDA)
  return cudaDeviceSetLimit(cudaLimitStackSize, bytes) == cudaSuccess;
#else
  (void)bytes;
  return false;
#endif
}
