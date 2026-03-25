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

/*! \file lapack_calls.cpp
 *
 * \brief Implementation of the interface with LAPACK libraries.
 */
#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

/*
#ifdef USE_LAPACK
#ifdef USE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
*/

#ifdef USE_LAPACK

#ifndef INCLUDE_LAPACK_CALLS_H_
#include "../include/lapack_calls.h"
#endif

#include <limits>

#ifdef USE_MKL
  extern "C" void zcopy_(np_int *n, MKL_Complex16 *arr1, np_int *inc1, MKL_Complex16 *arr2,
		     np_int *inc2);
  extern "C" void zgemm_(char *transa, char *transb, np_int *l, np_int *m, np_int *n,
		     MKL_Complex16 *alpha, MKL_Complex16 *a, np_int *lda,
		     MKL_Complex16 *b, np_int *ldb, MKL_Complex16 *beta,
		     MKL_Complex16 *c, np_int *ldc);
  extern "C" void zaxpy_(np_int *n, MKL_Complex16 *alpha, MKL_Complex16 *arr1, np_int *inc1,
		     MKL_Complex16 *arr2, np_int *inc2);
  extern "C" np_int izamax_(np_int *n, MKL_Complex16 *arr1, np_int *inc1);
#else
  extern "C" void zcopy_(np_int *n, dcomplex *arr1, np_int *inc1, dcomplex *arr2, np_int *inc2);
  extern "C" void zgemm_(char *transa, char *transb, np_int *l, np_int *m, np_int *n,
		     dcomplex *alpha, dcomplex *a, np_int *lda,
		     dcomplex *b, np_int *ldb, dcomplex *beta,
		     dcomplex *c, np_int *ldc);
  extern "C" void zaxpy_(np_int *n, dcomplex *alpha, dcomplex *arr1, np_int *inc1,
		     dcomplex *arr2, np_int *inc2);
  extern "C" np_int izamax_(np_int *n, dcomplex *arr1, np_int *inc1);
#endif

void zinvert(dcomplex **mat, np_int n, int &jer) {
  jer = 0;
  dcomplex *arr = &(mat[0][0]);
  const dcomplex uim = 0.0 + 1.0 * I;

#ifdef USE_MKL
  MKL_Complex16 *arr2 = (MKL_Complex16 *) arr;
#endif
  
  np_int* IPIV = new np_int[n]();
  
#ifdef USE_MKL
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, arr2, n, IPIV);
  LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, arr2, n, IPIV);
#else
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, arr, n, IPIV);
  LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, arr, n, IPIV);
#endif

  delete[] IPIV;
}

void zinvert_and_refine(dcomplex **mat, np_int n, int &jer, int &maxiters, double &accuracygoal, int refinemode) {

  jer = 0;
#ifdef USE_MKL
  MKL_Complex16 *arr = (MKL_Complex16 *) &(mat[0][0]);
#else
  dcomplex *arr = &(mat[0][0]);
#endif
  np_int nn = n*n;
  np_int incx = 1;
  np_int incx0 = 0;
#ifdef USE_MKL
  MKL_Complex16 *arr_orig = NULL;
  MKL_Complex16 *arr_residual = NULL;
  MKL_Complex16 *arr_refine = NULL;
  MKL_Complex16 *id = NULL;
#else
  dcomplex *arr_orig = NULL;
  dcomplex *arr_residual = NULL;
  dcomplex *arr_refine = NULL;
  dcomplex *id = NULL;
#endif
  if (maxiters>0) {
#ifdef USE_MKL
    arr_orig = new MKL_Complex16[nn];
    arr_residual = new MKL_Complex16[nn];
    arr_refine = new MKL_Complex16[nn];
    id = new MKL_Complex16[1];
    id[0].real =  1;
    id[0].imag =  0;
#else
    arr_orig = new dcomplex[nn];
    arr_residual = new dcomplex[nn];
    arr_refine = new dcomplex[nn];
    id = new dcomplex[1];
    id[0] = (dcomplex) 1;
#endif
    zcopy_(&nn, arr, &incx, arr_orig, &incx);
  }
  // const dcomplex uim = 0.0 + 1.0 * I;
  
  np_int* IPIV = new np_int[n]();
  
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, arr, n, IPIV);
  LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, arr, n, IPIV);
  delete[] IPIV;

  if (maxiters>0) {
    bool iteraterefine = true;
    char transa = 'N';
#ifdef USE_MKL
    MKL_Complex16 dczero;
    dczero.real = 0;
    dczero.imag = 0;
    MKL_Complex16 dcone;
    dcone.real = 1;
    dcone.imag = 0;
    MKL_Complex16 dcmone;
    dcmone.real = -1;
    dcmone.imag = 0;
#else
    dcomplex dczero = 0;
    dcomplex dcone = 1;
    dcomplex dcmone = -1;
#endif
    // multiply minus the original matrix times the inverse matrix
    // NOTE: factors in zgemm are swapped because zgemm is designed for column-major
    // Fortran-style arrays, whereas our arrays are C-style row-major.
    zgemm_(&transa, &transa, &n, &n, &n, &dcmone, arr, &n, arr_orig, &n, &dczero, arr_residual, &n);
    np_int incy = n+1;
    zaxpy_(&n, &dcone, id, &incx0, arr_residual, &incy);
    double oldmax = 0;
    if (refinemode >0) {
      np_int maxindex = izamax_(&nn, arr_residual, &incx);
#ifdef USE_MKL
      oldmax = cabs(arr_residual[maxindex].real + I*arr_residual[maxindex].imag);
#else
      oldmax = cabs(arr_residual[maxindex]);
#endif
      printf("Initial max residue = %g\n", oldmax);
      if (oldmax < accuracygoal) iteraterefine = false;
    }
    int iter;
    for (iter=0; (iter<maxiters) && iteraterefine; iter++) {
      zgemm_(&transa, &transa, &n, &n, &n, &dcone, arr_residual, &n, arr, &n, &dczero, arr_refine, &n);
      zaxpy_(&nn, &dcone, arr_refine, &incx, arr, &incx);
	// zcopy_(&nn, arr_refine, &incx, arr, &incx);
      zgemm_(&transa, &transa, &n, &n, &n, &dcmone, arr, &n, arr_orig, &n, &dczero, arr_residual, &n);
      zaxpy_(&n, &dcone, id, &incx0, arr_residual, &incy);
      if ((refinemode==2) || ((refinemode==1) && (iter == (maxiters-1)))) {
	np_int maxindex = izamax_(&nn, arr_residual, &incx);
#ifdef USE_MKL
	double newmax = cabs(arr_residual[maxindex].real + I*arr_residual[maxindex].imag);
#else
	double newmax = cabs(arr_residual[maxindex]);
#endif
	printf("Max residue after %d iterations = %g\n", iter+1, newmax);
	if ((refinemode==2) && ((newmax > oldmax)||(newmax < accuracygoal))) iteraterefine = false;
	oldmax = newmax; 
      }
    }
    if (refinemode==2) maxiters = iter;
    accuracygoal = oldmax;
    delete[] id;
    delete[] arr_refine;
    delete[] arr_orig;
    delete[] arr_residual;
  }

}

#endif
