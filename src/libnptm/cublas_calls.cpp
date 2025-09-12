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

/*! \file cublas_calls.cpp
 *
 * \brief Implementation of the interface with CUBLAS libraries.
 */
#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

//#define USE_CUBLAS 1
#ifdef USE_CUBLAS

#ifndef INCLUDE_CUBLAS_CALLS_H_
#include "../include/cublas_calls.h"
#endif

#include <limits>
#include <cuda_runtime.h>
#include <cublas_v2.h>

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

void cublas_zinvert(dcomplex **mat, np_int n, int device_id) {
  cudacall(cudaSetDevice(device_id));
  cublasHandle_t handle;
  cublascall(cublasCreate_v2(&handle));
  int batchsize = 1;
  int *piv, *info; // array of pivot indices
  np_int m = (np_int) n; // changed rows; a - mxm matrix
  np_int mm = m * m; // size of a
  cudacall(cudaMalloc<int>(&piv, m*batchsize*sizeof(int)));
  cudacall(cudaMalloc<int>(&info, batchsize*sizeof(int)));
  cuDoubleComplex *a = (cuDoubleComplex *)&(mat[0][0]); // pointer to first element on host
  cuDoubleComplex *d_a; // pointer to first element on device
  cudacall(cudaMalloc<cuDoubleComplex>(&d_a,m*m*sizeof(cuDoubleComplex)));
  cudacall(cudaMemcpy(d_a, a, m*m*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice));
  cuDoubleComplex **batch_d_a;
  cudacall(cudaMalloc<cuDoubleComplex*>(&batch_d_a,batchsize*sizeof(cuDoubleComplex*)));
  cudacall(cudaMemcpy(batch_d_a, &d_a, batchsize*sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice));
  cublascall(cublasZgetrfBatched(handle, m, batch_d_a, m, piv, info, batchsize));
  cuDoubleComplex *d_c; // this will contain the inverted matrix on the device
  cudacall(cudaMalloc<cuDoubleComplex>(&d_c,m*m*sizeof(cuDoubleComplex)));
  cuDoubleComplex **batch_d_c;
  cudacall(cudaMalloc<cuDoubleComplex*>(&batch_d_c,batchsize*sizeof(cuDoubleComplex*)));
  cudacall(cudaMemcpy(batch_d_c, &d_c, batchsize*sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice));
  cublascall(cublasZgetriBatched(handle,n,batch_d_a,m,piv,batch_d_c,m,info,batchsize));
  cudacall(cudaMemcpy(a,d_c,m*m*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
  cudaFree(batch_d_a);
  cudaFree(batch_d_c);
  cudaFree(piv);
  cudaFree(info);
  cudaFree(d_a);
  cudaFree(d_c);
  cublasDestroy_v2(handle);

}

void cublas_zinvert_and_refine(dcomplex **mat, np_int n, int &maxiters, double &accuracygoal, int refinemode, int device_id) {
  cudacall(cudaSetDevice(device_id));
  cublasHandle_t handle;
  cublascall(cublasCreate_v2(&handle));
  int batchsize = 1;
  int *piv, *info; // array of pivot indices
  np_int m = (np_int) n; // changed rows; a - mxm matrix
  np_int mm = m * m; // size of a
  cudacall(cudaMalloc<int>(&piv, m*batchsize*sizeof(int)));
  cudacall(cudaMalloc<int>(&info, batchsize*sizeof(int)));
  cuDoubleComplex *a = (cuDoubleComplex *)&(mat[0][0]); // pointer to first element on host
  cuDoubleComplex *d_a; // pointer to first element on device
  cudacall(cudaMalloc<cuDoubleComplex>(&d_a,m*m*sizeof(cuDoubleComplex)));
  cudacall(cudaMemcpy(d_a, a, m*m*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice));
  cuDoubleComplex **batch_d_a;
  cudacall(cudaMalloc<cuDoubleComplex*>(&batch_d_a,batchsize*sizeof(cuDoubleComplex*)));
  cudacall(cudaMemcpy(batch_d_a, &d_a, batchsize*sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice));
  cublascall(cublasZgetrfBatched(handle, m, batch_d_a, m, piv, info, batchsize));
  cuDoubleComplex *d_c; // this will contain the inverted matrix on the device
  cudacall(cudaMalloc<cuDoubleComplex>(&d_c,m*m*sizeof(cuDoubleComplex)));
  cuDoubleComplex **batch_d_c;
  cudacall(cudaMalloc<cuDoubleComplex*>(&batch_d_c,batchsize*sizeof(cuDoubleComplex*)));
  cudacall(cudaMemcpy(batch_d_c, &d_c, batchsize*sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice));
  cublascall(cublasZgetriBatched(handle,m,batch_d_a,m,piv,batch_d_c,m,info,batchsize));
  //cudacall(cudaMemcpy(a,d_c,m*m*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
  cudaFree(batch_d_a);
  cudaFree(batch_d_c);
  cudaFree(piv);
  cudaFree(info);
  cuDoubleComplex *d_a_residual;
  cuDoubleComplex *d_a_refine;
  cuDoubleComplex *d_id;
  if (maxiters>0) {
    // copy the original matrix again to d_a, so I do not need to destroy the old d_a and recreate a new one
    cudacall(cudaMemcpy(d_a, a, m*m*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice)); // from here on, d_a contains the original matrix, for refinement use  
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
  }
  bool iteraterefine = true;
  if (maxiters>0) {
    cuDoubleComplex cu_mone;
    (((double *) &(cu_mone))[0]) = -1;
    (((double *) &(cu_mone))[1]) = 0;
    cuDoubleComplex cu_one;
    (((double *) &(cu_one))[0]) = 1;
    (((double *) &(cu_one))[1]) = 0;
    cuDoubleComplex cu_zero;
    (((double *) &(cu_zero))[0]) = 0;
    (((double *) &(cu_zero))[1]) = 0;
    // multiply minus the original matrix times the inverse matrix
    // NOTE: factors in zgemm are swapped because zgemm is designed for column-major
    // Fortran-style arrays, whereas our arrays are C-style row-major.
    cublascall(CUZGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_mone, d_c, m, d_a, m, &cu_zero, d_a_residual, m));
    // add the identity to the product
    cublascall(CUZAXPY(handle, m, &cu_one, d_id, 0, d_a_residual, m+1));
    double oldmax=0;
    if (refinemode >0) {
      np_int maxindex;
      // find the maximum absolute value of the residual
      cublascall(CUIZAMAX(handle, mm, d_a_residual, 1, &maxindex));
      cuDoubleComplex cublasmax;
      // transfer the maximum value to the host
      cudacall(cudaMemcpy(&cublasmax, d_a_residual+maxindex, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
      // take the module
      oldmax = cabs( (((double *) &(cublasmax))[0]) + I*(((double *) &(cublasmax))[1]));
      printf("Initial max residue = %g\n", oldmax);
      if (oldmax < accuracygoal) iteraterefine = false;
    }
    // begin correction loop (should iterate maxiters times)
    int iter;
    for (iter=0; (iter<maxiters) && iteraterefine; iter++) {
      // multiply the inverse times the residual, add to the initial inverse
      cublascall(CUZGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_one, d_a_residual, m, d_c, m, &cu_zero, d_a_refine, m));
      // add to the initial inverse
      cublascall(CUZAXPY(handle, mm, &cu_one, d_a_refine, 1, d_c, 1));
      // multiply minus the original matrix times the new inverse matrix
      cublascall(CUZGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, m, m, &cu_mone, d_c, m, d_a, m, &cu_zero, d_a_residual, m));
      // add the identity to the product
      cublascall(CUZAXPY(handle, m, &cu_one, d_id, 0, d_a_residual, m+1));
      if ((refinemode==2) || ((refinemode==1) && (iter == (maxiters-1)))) {
	// find the maximum absolute value of the residual
	np_int maxindex;
	cublascall(CUIZAMAX(handle, mm, d_a_residual, 1, &maxindex));
	// transfer the maximum value to the host
	cuDoubleComplex cublasmax;
	cudacall(cudaMemcpy(&cublasmax, d_a_residual+maxindex, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
	// take the module
	double newmax = cabs( (((double *) &(cublasmax))[0]) + I*(((double *) &(cublasmax))[1]));
	printf("Max residue after %d iterations = %g\n", iter+1, newmax);
	// if the maximum in the residual decreased from the previous iteration,
	// update oldmax and go on, otherwise no point further iterating refinements
	if ((refinemode==2) && ((newmax > oldmax)||(newmax < accuracygoal))) iteraterefine = false;
	oldmax = newmax;
      }
    }
    // if we are being called with refinemode=2, then on exit we set maxiters to the actual number of iters we performed to achieve the required accuracy
    if (refinemode==2) maxiters = iter;
    accuracygoal = oldmax;
    // end correction loop
  }
  // copy the final refined inverted matrix back from device to host
  cudacall(cudaMemcpy(a,d_c,m*m*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost));
  // free temporary device arrays 
  cudaFree(d_a);
  cudaFree(d_c);
  if (maxiters>0) {
    cudaFree(d_id);
    cudaFree(d_a_refine);
    cudaFree(d_a_residual);
  }

  cublasDestroy_v2(handle);

}


#endif
