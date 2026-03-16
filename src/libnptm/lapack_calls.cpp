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


#ifdef USE_LAPACK

#include <string>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifdef USE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif


#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_LAPACK_CALLS_H_
#include "../include/lapack_calls.h"
#endif

extern "C" void zcopy_(np_int *n, lcomplex *arr1, np_int *inc1, lcomplex *arr2, np_int *inc2);
extern "C" void zgemm_(
  char *transa, char *transb, np_int *l, np_int *m, np_int *n, lcomplex *alpha, lcomplex *a,
  np_int *lda, lcomplex *b, np_int *ldb, lcomplex *beta, lcomplex *c, np_int *ldc
);
extern "C" void zaxpy_(
  np_int *n, lcomplex *alpha, lcomplex *arr1, np_int *inc1, lcomplex *arr2, np_int *inc2
);
extern "C" np_int izamax_(np_int *n, lcomplex *arr1, np_int *inc1);

using namespace std;

void zinvert(dcomplex **mat, np_int n, int &jer, const RuntimeSettings& rs) {
  jer = 0;
  char buffer[128];
  string message;
  lapack_int info, inc1 = 1;
  lcomplex *arr = &(mat[0][0]);
  lcomplex *arr_orig;
  lcomplex lapack_one = 1.0 + I * 0.0;
  np_int nn = n * n;
  if (rs.use_refinement && rs.invert_mode != RuntimeSettings::INV_MODE_RBT) {
    lapack_int inc1 = 1;
    arr_orig = new lcomplex[nn];
    zcopy_(&nn, arr, &inc1, arr_orig, &inc1);
  }
  if (rs.invert_mode == RuntimeSettings::INV_MODE_LU) {
    // >>> LU INVERSION SECTION <<<
    np_int* IPIV = new np_int[n];
    LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, arr, n, IPIV);
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, arr, n, IPIV);
    delete[] IPIV;
    if (rs.use_refinement) {
      info = lapack_newton(rs, arr_orig, n, arr);
    }
    // >>> END OF LU INVERSION SECTION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_GESV) {
    // >>> GESV INVERSION SECTION <<<
    lcomplex *id = new lcomplex[nn]();
    lapack_int *piv = new lapack_int[n];
    for (lapack_int i = 0; i < n; i++) {
      id[i * (n + 1)] = lapack_one;
    }
    LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, n, arr, n, piv, id, n);
    if (info != LAPACK_SUCCESS) {
      message = "ERROR: call to zgesv_() returned info code " + to_string(info) + "!\n";
      rs.logger->err(message);
      exit(1);
    }
    delete[] piv; // free host memory
    if (rs.use_refinement) {
      info = lapack_newton(rs, arr_orig, n, id);
    }
    zcopy_(&nn, id, &inc1, arr, &inc1);
    delete[] id;
    // >>> END OF GESV INVERSION SECTION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_RBT) {
    // >>> RBT INVERSION SECTION <<<
    // RBT inversion not implemented
    message = "ERROR: not implemented!\n";
    rs.logger->err(message);
    exit(1);
    // >>> END OF RBT INVERSION SECTION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_SVD) {
    // >>> SVD INVERSION SECTION <<<
    // SVD inversion not implemented
    message = "ERROR: not implemented!\n";
    rs.logger->err(message);
    exit(1);
    // >>> END OF SVD INVERSION SECTION <<<
  } // inversion mode switch
  if (rs.use_refinement && rs.invert_mode != RuntimeSettings::INV_MODE_RBT) {
    delete[] arr_orig;
  }
}

lapack_int lapack_newton(
  const RuntimeSettings& rs, lcomplex* a_orig, lapack_int m, lcomplex* a
) {
  lapack_int err = LAPACK_SUCCESS;
  string message;
  char buffer[128];
  char lapackNoTrans = 'N';
  const int max_ref_iters = rs.max_ref_iters;
  lcomplex lapack_zero = 0.0 + I * 0.0;
  lcomplex lapack_one = 1.0 + I * 0.0;
  lcomplex lapack_mone = -1.0 + I * 0.0;
  lapack_int mm = m * m;
  lapack_int incx, incy;
  lcomplex *ax, *r;
  lcomplex *id_diag = new lcomplex[m];
  double oldmax = 2.0e+16, curmax = 1.0e+16;
  for (lapack_int hi = 0; hi < m; hi++)
    id_diag[hi] = lapack_one;
  ax = new lcomplex[mm];
  r = new lcomplex[mm];
  double max_residue, target_residue;
  incx = 1;
  lapack_int maxindex = izamax_(&mm, a, &incx) - 1;
  lcomplex lapackmax = a[maxindex];
  curmax = cabs(lapackmax); //cabs(magmamax.x + I * magmamax.y);
  target_residue = curmax * rs.accuracy_goal;
  sprintf(buffer, "INFO: largest matrix value has modulus %.5le; target residue is %.5le.\n", curmax, target_residue);
  message = buffer;
  rs.logger->log(message);
  for (int ri = 0; ri < max_ref_iters; ri++) {
    oldmax = curmax;
    // Compute -A*X
    zgemm_(
      &lapackNoTrans, &lapackNoTrans, &m, &m, &m, &lapack_mone, a, &m,
      a_orig, &m, &lapack_zero, ax, &m
    );
    // Transform -A*X into (I - A*X)
    incx = 1;
    incy = m + 1;
    zaxpy_(&m, &lapack_one, id_diag, &incx, ax, &incy);
    maxindex = izamax_(&mm, ax, &incx) - 1;
    lapackmax = ax[maxindex];
    curmax = cabs(lapackmax);
    sprintf(buffer, "DEBUG: iteration %d has residue %.5le; target residue is %.5le.\n", ri, curmax, target_residue);
    message = buffer;
    rs.logger->log(message, LOG_DEBG);
    if (curmax < 0.99 * oldmax) {
      // Compute R = (I - A*X)*X
      zgemm_(
        &lapackNoTrans, &lapackNoTrans, &m, &m, &m, &lapack_one, a, &m,
	ax, &m, &lapack_zero, r, &m
      );
      // Set X = X + R
      zaxpy_(&mm, &lapack_one, r, &incx, a, &incx);
      if (curmax < target_residue) {
	message = "DEBUG: good news - optimal convergence achieved. Stopping.\n";
	rs.logger->log(message, LOG_DEBG);
	break; // ri for
      }
    } else {
      message = "WARN: not so good news - cannot improve further. Stopping.\n";
      rs.logger->log(message, LOG_WARN);
      break; // ri for
    }
  }
  delete[] id_diag;
  delete[] ax;
  delete[] r;
  return err;
}

#endif
