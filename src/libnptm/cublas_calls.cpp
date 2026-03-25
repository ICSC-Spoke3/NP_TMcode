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
 */

#ifdef USE_CUBLAS

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <limits>
#include <string>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
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
#define CUZGEMM cublasZgemm_64
#define CUZAXPY cublasZaxpy_64
#define CUIZAMAX cublasIzamax_64
#else
#define CUZGEMM cublasZgemm
#define CUZAXPY cublasZaxpy
#define CUIZAMAX cublasIzamax
#endif

#define cudacall(call)							\
  do									\
    {									\
      cudaError_t err = (call);						\
      if(cudaSuccess != err)						\
        {								\
	  fprintf(stderr,"CUDA Error:\nFile = %s\nLine = %d\nReason = %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
	  cudaDeviceReset();						\
	  exit(EXIT_FAILURE);						\
        }								\
    }									\
  while (0)

#define cublascall(call)						\
  do									\
    {									\
      cublasStatus_t status = (call);					\
      if(CUBLAS_STATUS_SUCCESS != status)				\
        {								\
	  fprintf(stderr,"CUBLAS Error:\nFile = %s\nLine = %d\nCode = %d\n", __FILE__, __LINE__, status); \
	  cudaDeviceReset();						\
	  exit(EXIT_FAILURE);						\
        }								\
      									\
    }									\
  while(0)

using namespace std;

void cublas_zinvert(dcomplex **mat, np_int n, int device_id, const RuntimeSettings& rs) {
  char buffer[128];
  string message;
  int ref_err = 0;
  cudacall(cudaSetDevice(device_id));
  cublasHandle_t handle;
  cublascall(cublasCreate_v2(&handle));
  int batchsize = 1;
  if (rs.invert_mode == RuntimeSettings::INV_MODE_LU) {
    int *piv, *info; // array of pivot indices
    np_int m = (np_int) n; // changed rows; a - mxm matrix
    np_int mm = m * m; // size of a
    cudacall(cudaMalloc<int>(&piv, m*batchsize*sizeof(int)));
    cudacall(cudaMalloc<int>(&info, batchsize*sizeof(int)));
    cuDoubleComplex *a = (cuDoubleComplex *)&(mat[0][0]); // pointer to first element on host
    cuDoubleComplex *d_a; // pointer to first element on device
    cudacall(cudaMalloc<cuDoubleComplex>(&d_a,m*m*sizeof(cuDoubleComplex)));
    cudacall(cudaMemcpy(d_a, a, m*m*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice)); // copy a -> d_a
    cuDoubleComplex **batch_d_a;
    cudacall(cudaMalloc<cuDoubleComplex*>(&batch_d_a,batchsize*sizeof(cuDoubleComplex*)));
    cudacall(cudaMemcpy(batch_d_a, &d_a, batchsize*sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice));
    cublascall(cublasZgetrfBatched(handle, m, batch_d_a, m, piv, info, batchsize)); // call to ZGETRF
    cuDoubleComplex *d_c; // this will contain the inverted matrix on the device
    cudacall(cudaMalloc<cuDoubleComplex>(&d_c,m*m*sizeof(cuDoubleComplex)));
    cuDoubleComplex **batch_d_c;
    cudacall(cudaMalloc<cuDoubleComplex*>(&batch_d_c,batchsize*sizeof(cuDoubleComplex*)));
    cudacall(cudaMemcpy(batch_d_c, &d_c, batchsize*sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice));
    cublascall(cublasZgetriBatched(handle,n,batch_d_a,m,piv,batch_d_c,m,info,batchsize));
    if (rs.use_refinement) {
      cuDoubleComplex cu_mone;
      (((double *) &(cu_mone))[0]) = -1;
      (((double *) &(cu_mone))[1]) = 0;
      cuDoubleComplex cu_one;
      (((double *) &(cu_one))[0]) = 1;
      (((double *) &(cu_one))[1]) = 0;
      cuDoubleComplex cu_zero;
      (((double *) &(cu_zero))[0]) = 0;
      (((double *) &(cu_zero))[1]) = 0;

      cuDoubleComplex *d_a_residual;
      cuDoubleComplex *d_a_refine;
      cuDoubleComplex *d_id;
      // copy the original matrix again to d_a, so I do not need to destroy the old d_a and recreate a new one
      cudacall(cudaMemcpy(d_a, a, m*m*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice)); // from here on, d_a contains the original matrix, for refinement use  
      // copy the unrefined inverted matrix from device to host, so we have a back-up if refinement fails
      cudacall(cudaMemcpy(a, d_c, m*m*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
      cudacall(cudaMalloc<cuDoubleComplex>(&d_a_residual, m*m*sizeof(cuDoubleComplex)));
      cudacall(cudaMalloc<cuDoubleComplex>(&d_a_refine, m*m*sizeof(cuDoubleComplex)));
      // allocate memory for the temporary matrix products
      dcomplex *native_id = new dcomplex[1];
      native_id[0] = 1;
      cuDoubleComplex *m_id = (cuDoubleComplex *) &(native_id[0]);
      // fill it with 1
      cudacall(cudaMalloc<cuDoubleComplex>(&d_id, sizeof(cuDoubleComplex)));
      cudacall(cudaMemcpy(d_id, m_id, sizeof(cuDoubleComplex),cudaMemcpyHostToDevice)); // copy identity to device vector
      delete[] native_id; // free identity vector on host

      // Detect the maximum value of the inverse matrix.
      double oldmax = 2.0e+16, curmax = 1.0e+16;
      np_int maxindex;
      cuDoubleComplex cublasmax;
      cublascall(CUIZAMAX(handle, mm, d_c, 1, &maxindex));
      cudacall(cudaMemcpy(&cublasmax, d_c + maxindex - 1, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
      curmax = cabs( (((double *) &(cublasmax))[0]) + I * (((double *) &(cublasmax))[1]));
      sprintf(buffer, "INFO: largest matrix value has modulus %.5le.\n", curmax);
      message = buffer;
      rs.logger->log(message);

      // Iterative refinement loop
      int max_ref_iters = rs.max_ref_iters;
      for (int ri = 0; ri < max_ref_iters; ri++) {
	oldmax = curmax;
	// multiply minus the original matrix times the inverse matrix
	// NOTE: factors in zgemm are swapped because zgemm is designed for column-major
	// Fortran-style arrays, whereas our arrays are C-style row-major.
	cublascall(CUZGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_mone, d_c, m, d_a, m, &cu_zero, d_a_residual, m));
	// add the identity to the product
	cublascall(CUZAXPY(handle, m, &cu_one, d_id, 0, d_a_residual, m+1));
	// find the maximum absolute value of the residual
	cublascall(CUIZAMAX(handle, mm, d_a_residual, 1, &maxindex));
	// transfer the maximum value to the host
	cudacall(cudaMemcpy(&cublasmax, d_a_residual + maxindex - 1, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
	// take the module
	curmax = cabs( (((double *) &(cublasmax))[0]) + I * (((double *) &(cublasmax))[1]));
	sprintf(buffer, "DEBUG: iteration %d has residue %.5le; target residue is %.5le.\n", ri, curmax, rs.accuracy_goal);
	message = buffer;
	rs.logger->log(message, LOG_DEBG);
	if (curmax < 0.5) { // Safe conditions for Newton-Schultz iteration.
	  if (curmax < 0.95 * oldmax) { // Newton-Schultz iteration is improving and can proceed.
	    // multiply the inverse times the residual, add to the initial inverse
	    cublascall(CUZGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_one, d_a_residual, m, d_c, m, &cu_zero, d_a_refine, m));
	    // add to the initial inverse
	    cublascall(CUZAXPY(handle, mm, &cu_one, d_a_refine, 1, d_c, 1));
	    if (curmax < rs.accuracy_goal) {
	      message = "DEBUG: good news - optimal convergence achieved. Stopping.\n";
	      rs.logger->log(message, LOG_DEBG);
	      ref_err = 0;
	      break; // ri for
	    }
	  } else {
	    if (curmax > 0.1) {
	      sprintf(buffer, "INFO: iteration %d achieved limiting residue %.5le. Cannot reach goal. Reverting.\n", ri, curmax);
	      message = buffer;
	      rs.logger->log(message);
	      ref_err = 1;
	    } else {
	      sprintf(buffer, "WARN: iteration %d achieved limiting residue %.5le. Stopping.\n", ri, curmax);
	      message = buffer;
	      rs.logger->log(message, LOG_WARN);
	    }
	    break; // ri for
	  }
	} else { // curmax > 0.5. Newton-Schultz iteration is dangerous.
	  if (curmax < oldmax) {
	    sprintf(buffer, "WARN: iteration %d has residue %.5le. Iterating is dangerous. Stopping.\n", ri, curmax);
	    message = buffer;
	    rs.logger->log(message, LOG_WARN);
	  } else {
	    ref_err = 1;
	    sprintf(buffer, "INFO: iteration %d has residue %.5le. Reverting to unrefined and stopping.\n", ri, curmax);
	    message = buffer;
	    rs.logger->log(message);
	  }
	  break; // ri for
	}
      } // end of ri for
      cudaFree(d_a_residual);
      cudaFree(d_a_refine);
      cudaFree(d_id);
    }
    if (ref_err == 0) { // Refinement was fine, we can copy refined matrix back to host.
      cudacall(cudaMemcpy(a, d_c, m * m * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    }
    cudaFree(batch_d_a);
    cudaFree(batch_d_c);
    cudaFree(piv);
    cudaFree(info);
    cudaFree(d_a);
    cudaFree(d_c);
  } else {
    message = "ERROR: cuBLAS solver only implements LU inversion!\n";
    rs.logger->err(message);
    exit(1);
  }
  cublasDestroy_v2(handle);
}

#endif
