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

/*! \file magma_calls.cpp
 *
 * \brief Implementation of the interface with MAGMA libraries.
 */

#ifdef USE_MAGMA

//#include <limits>
#include <cstdio>
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

#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif

#ifdef DEBUG_REFINE
#include <hdf5.h>
#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif // INCLUDE_ERRORS_H_
#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif // INCLUDE_LIST_H_
#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif // INCLUDE_FILE_IO_H_
#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif // INCLUDE_UTILS_H_

#endif // DEBUG_REFINE

using namespace std;

void magma_zinvert(dcomplex **mat, np_int n, int &jer, int device_id, const RuntimeSettings& rs) {
  magmaDoubleComplex *a = (magmaDoubleComplex *)mat[0]; // pointer to first element on host
  string message;
  char buffer[128];
  magma_int_t err = MAGMA_SUCCESS;
  magma_queue_t queue = NULL;
  magma_device_t dev = (magma_device_t)device_id;
  magma_queue_create(dev, &queue);
  magma_int_t info;
  magma_int_t m = (magma_int_t)n; // changed rows; a - mxm matrix
  magma_int_t mm = m * m; // size of a
  const magmaDoubleComplex magma_zero = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex magma_one = MAGMA_Z_MAKE(1.0, 0.0);
  const magmaDoubleComplex magma_mone = MAGMA_Z_MAKE(-1.0, 0.0);
  const magmaDoubleComplex magma_two = MAGMA_Z_MAKE(2.0, 0.0);
  if (rs.invert_mode == RuntimeSettings::INV_MODE_LU) { 
    // >>> LU INVERSION <<<
    magmaDoubleComplex *d_a; // pointer to first element on device
    magmaDoubleComplex *dwork; // work space pointer on device
    magma_int_t ldwork = m * magma_get_zgetri_nb(m); // optimal block size
    err = magma_zmalloc(&dwork, ldwork); // dev. mem. for ldwork
    magma_int_t *piv = new magma_int_t[m];
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate dwork!\n";
      rs.logger->err(message);
      exit(1);
    }
    err = magma_zmalloc(&d_a, mm); // dev. mem. for ldwork
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
    magma_zgetrf_gpu(m, m, d_a, m, piv, &info);
    magma_zgetri_gpu(m, d_a, m, piv, dwork, ldwork, &info);
    delete[] piv; // delete piv created by magma_zgetrf()
    magma_free(dwork);
    if (rs.use_refinement) {
      err = magma_newton(rs, a, m, d_a, queue);
    }
    magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy d_a -> a
    magma_free(d_a);
    // >>> END OF LU INVERSION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_GESV) {
    // >>> GESV INVERSION <<<
    magmaDoubleComplex *h_id = new magmaDoubleComplex[mm]();
    magmaDoubleComplex *d_a, *d_id;
    magma_int_t *piv = new magma_int_t[m];
    for (magma_int_t i = 0; i < m; i++) {
      h_id[i * (m + 1)] = magma_one;
    }
    err = magma_zmalloc(&d_id, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_id!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, h_id, m, d_id , m, queue); // copy identity matrix to device
    delete[] h_id;
    // allocate matrix on device
    err = magma_zmalloc(&d_a, mm); // device memory for a
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
    magma_zgesv_gpu(m, m, d_a, m, piv, d_id, m, &info);
    delete[] piv; // free host memory
    magma_free(d_a);
    if (rs.use_refinement) {
      err = magma_newton(rs, a, m, d_id, queue);
    }
    magma_zgetmatrix(m, m, d_id , m, a, m, queue); // copy d_id -> a
    magma_free(d_id);
    // >>> END OF GESV INVERSION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_RBT) {
    // >>> RBT INVERSION <<<
    magmaDoubleComplex *d_a; // pointer to first element on device
    magma_bool_t iterate = (magma_bool_t)rs.use_refinement;
    magma_int_t *piv = new magma_int_t[m]; // host mem.
    magmaDoubleComplex *h_id = new magmaDoubleComplex[mm]();
    for (np_int i = 0; i < m; i++) {
      h_id[i * (m + 1)] = magma_one;
    }
    err = magma_zmalloc(&d_a, mm); // device memory for a
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
    magma_zgesv_rbt(iterate, m, m, a, m, h_id, m, &info);
    memcpy(a, h_id, mm * sizeof(magmaDoubleComplex));
    delete[] h_id;
    delete[] piv;
    magma_free(d_a); 
    // >>> END OF RBT INVERSION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_SVD) {
    // >>> SVD INVERSION <<<
    // Step 1: compute singular value decomposition
    magma_int_t lwork;
    magmaDoubleComplex *h_work = new magmaDoubleComplex[m];
    magma_zgesdd(
      MagmaAllVec, m, m, NULL, m, NULL, NULL, m, NULL, m, h_work, -1, NULL, NULL, &info
    );
    if (info == 0) {
      lwork = (magma_int_t)(MAGMA_Z_REAL(h_work[0]));
      message = "DEBUG: lwork = " + to_string(lwork) + "\n";
      rs.logger->log(message, LOG_DEBG);
    } else {
      message = "ERROR: magma_zgesdd failed on resource querying call!\n";
      rs.logger->err(message);
      exit(1);
    }
    delete[] h_work;
    h_work = new magmaDoubleComplex[lwork];
    magmaDoubleComplex *h_a, *h_U, *h_Vt;
    h_a = new magmaDoubleComplex[mm];
    memcpy(h_a, a, mm * sizeof(magmaDoubleComplex));
    h_U = new magmaDoubleComplex[mm];
    h_Vt = new magmaDoubleComplex[mm];
    double *h_s = new double[m];
    double *h_rwork = new double[5 * mm + 7 * m];
    magma_int_t *h_iwork = new magma_int_t[8 * m];
    magma_zgesdd(
      MagmaAllVec, m, m, h_a, m, h_s, h_U, m, h_Vt, m, h_work, lwork, h_rwork, h_iwork, &info
    );
    delete[] h_work;
    delete[] h_rwork;
    delete[] h_iwork;
    delete[] h_a;
    sprintf(buffer, "%.5le", h_s[0]);
    message = "DEBUG: s[0] = ";
    message += buffer;
    message += "; s[";
    message += to_string(m - 1);
    message += "] = ";
    sprintf(buffer, "%.5le", h_s[m - 1]);
    message += buffer;
    message += ".\n";
    rs.logger->log(message, LOG_DEBG);
    
    // Step 2: Upload decomposed matix to GPU
    magmaDoubleComplex *d_U, *d_Vt, *d_a;
    err = magma_zmalloc(&d_Vt, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_Vt!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, h_Vt, m, d_Vt, m, queue);
    err = magma_zmalloc(&d_U, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_U!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, h_U, m, d_U, m, queue);
    err = magma_zmalloc(&d_a, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }

    // Step 3: compute threshold-truncated pseudo-inverse on GPU
    double eps_mach = 2.22e-18;
    double threshold = h_s[0] * m * eps_mach; 
    for (magma_int_t i = 0; i < m; i++) {
      magmaDoubleComplex inv_s = (h_s[i] > threshold) ? 
	MAGMA_Z_MAKE(1.0 / h_s[i], 0.0) : 
	MAGMA_Z_MAKE(0.0, 0.0);
      magma_zscal(m, inv_s, d_Vt + i, m, queue);
    }
    delete[] h_s;
    // d_a = d_vt^H * d_u^H
    magmablas_zgemm(
      MagmaConjTrans, MagmaConjTrans, m, m, m, magma_one, d_Vt, m, d_U, m, 
      magma_zero, d_a, m, queue
    );

    // Step 4: release matrix decomposition
    magma_free(d_U);
    magma_free(d_Vt);
    delete[] h_U;
    delete[] h_Vt;

    // Step 5: refine inversion
    if (rs.use_refinement) {
      err = magma_newton(rs, a, m, d_a, queue);
    }

    // Step 6: get result back to host
    magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy d_a -> a
    magma_free(d_a);
    // >>> END OF SVD INVERSION <<<
  }
  magma_queue_destroy(queue); // destroy queue
  jer = (int)err;
}

magma_int_t magma_newton_norm(
  const RuntimeSettings& rs, magmaDoubleComplex* a, const magma_int_t m,
  magmaDoubleComplex* d_a, magma_queue_t queue
) {
  magma_int_t err;
  string message;
  char buffer[128];
  const int max_ref_iters = rs.max_ref_iters;
  const magmaDoubleComplex magma_zero = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex magma_one = MAGMA_Z_MAKE(1.0, 0.0);
  const magmaDoubleComplex magma_mone = MAGMA_Z_MAKE(-1.0, 0.0);
  const magmaDoubleComplex magma_two = MAGMA_Z_MAKE(2.0, 0.0);
  const magma_int_t mm = m * m;
  magmaDoubleComplex *d_a_orig, *d_ax, *d_axa;
  magmaDoubleComplex *h_id_diag = new magmaDoubleComplex[m];
  magmaDoubleComplex *d_id_diag;
  for (magma_int_t hi = 0; hi < m; hi++)
    h_id_diag[hi] = magma_one;
  err = magma_zmalloc(&d_id_diag, m);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_id_diag!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_zsetvector(m, h_id_diag, 1, d_id_diag, 1, queue);
  delete[] h_id_diag;
  double norm = 1.0e+16, old_norm = 2.0e+16;
  err = magma_zmalloc(&d_a_orig, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_a_orig!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_zsetmatrix(m, m, a, m, d_a_orig, m, queue); // copy a -> d_a_orig
  err = magma_zmalloc(&d_ax, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_ax!\n";
    rs.logger->err(message);
    exit(1);
  }
  err = magma_zmalloc(&d_axa, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_axa!\n";
    rs.logger->err(message);
    exit(1);
  }
  double normA = magma_dznrm2(mm, d_a_orig, 1, queue);
  if (normA == 0.0) normA = 1.0;
  for (int ri = 0; ri < max_ref_iters; ri++) {
    // Compute A*X
    magmablas_zgemm(
      MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_a_orig, m,
      d_a, m, magma_zero, d_ax, m, queue
    );
    // Compute A*X*A
    magmablas_zgemm(
      MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_ax, m,
      d_a_orig, m, magma_zero, d_axa, m, queue
    );
    // Now we use d_axa in place of R: d_axa = A*X*A - A
    magma_zaxpy(mm, magma_mone, d_a_orig, 1, d_axa, 1, queue);
    // Compute the norm of R
    norm = magma_dznrm2(mm, d_axa, 1, queue);
    double relative_res = norm / normA;
    sprintf(buffer, "%.5le", relative_res);
    message = "INFO: refinement iteration " + to_string(ri) + " achieved relative residual of "
      + buffer + ".\n";
    rs.logger->log(message);
    if (norm < 0.99 * old_norm) {
      // Transform d_ax in 2*I - A*X
      magma_zscal(mm, magma_mone, d_ax, 1, queue);
      magma_zaxpy(m, magma_two, d_id_diag, 1, d_ax, m + 1, queue);
      // Replace d_axa with X*(2*I - A*X)
      magma_zgemm(
        MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_a, m,
	d_ax, m, magma_zero, d_axa, m, queue
      );
      // Set d_a = X*(2*I - A*X)
      magma_zcopy(mm, d_axa, 1, d_a, 1, queue);
      old_norm = norm;
      if (relative_res < 1.0e-16) {
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
  magma_free(d_a_orig);
  magma_free(d_ax);
  magma_free(d_axa);
  magma_free(d_id_diag);
  return err;
}

magma_int_t magma_newton(
  const RuntimeSettings& rs, magmaDoubleComplex* a, const magma_int_t m,
  magmaDoubleComplex* d_a, magma_queue_t queue
) {
  magma_int_t err;
  string message;
  char buffer[128];
  const int max_ref_iters = rs.max_ref_iters;
  const magmaDoubleComplex magma_zero = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex magma_one = MAGMA_Z_MAKE(1.0, 0.0);
  const magmaDoubleComplex magma_mone = MAGMA_Z_MAKE(-1.0, 0.0);
  const magmaDoubleComplex magma_two = MAGMA_Z_MAKE(2.0, 0.0);
  const magma_int_t mm = m * m;
  magmaDoubleComplex *d_a_orig, *d_ax, *d_r;
  magmaDoubleComplex *h_id_diag = new magmaDoubleComplex[m];
  magmaDoubleComplex *d_id_diag;
  double oldmax = 2.0e+16, curmax = 1.0e+16;
  for (magma_int_t hi = 0; hi < m; hi++)
    h_id_diag[hi] = magma_one;
  err = magma_zmalloc(&d_id_diag, m);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_id_diag!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_zsetvector(m, h_id_diag, 1, d_id_diag, 1, queue);
  delete[] h_id_diag;
  err = magma_zmalloc(&d_a_orig, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_a_orig!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_zsetmatrix(m, m, a, m, d_a_orig, m, queue); // copy a -> d_a_orig
  err = magma_zmalloc(&d_ax, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_ax!\n";
    rs.logger->err(message);
    exit(1);
  }
  err = magma_zmalloc(&d_r, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_r!\n";
    rs.logger->err(message);
    exit(1);
  }
  double max_residue, target_residue;
  magma_int_t maxindex = magma_izamax(mm, d_a, 1, queue) - 1;
  magmaDoubleComplex magmamax = magma_mone;
  magma_zgetvector(1, d_a + maxindex, 1, &magmamax, 1, queue);
  curmax = MAGMA_Z_ABS(magmamax); //cabs(magmamax.x + I * magmamax.y);
  target_residue = curmax * rs.accuracy_goal;
  sprintf(buffer, "INFO: largest matrix value has modulus %.5le; target residue is %.5le.\n", curmax, target_residue);
  message = buffer;
  rs.logger->log(message);
  for (int ri = 0; ri < max_ref_iters; ri++) {
    oldmax = curmax;
    // Compute -A*X
    magmablas_zgemm(
      MagmaNoTrans, MagmaNoTrans, m, m, m, magma_mone, d_a_orig, m,
      d_a, m, magma_zero, d_ax, m, queue
    );
    // Transform -A*X into (I - A*X)
    magma_zaxpy(m, magma_one, d_id_diag, 1, d_ax, m + 1, queue);
    maxindex = magma_izamax(mm, d_ax, 1, queue) - 1;
    magma_zgetvector(1, d_ax + maxindex, 1, &magmamax, 1, queue);
    curmax = cabs(magmamax.x + I * magmamax.y);
    sprintf(buffer, "DEBUG: iteration %d has residue %.5le; target residue is %.5le.\n", ri, curmax, target_residue);
    message = buffer;
    rs.logger->log(message, LOG_DEBG);
    if (curmax < 0.99 * oldmax) {
      // Compute R = (A*X - I)*X
      magmablas_zgemm(
        MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_ax, m,
	d_a, m, magma_zero, d_r, m, queue
      );
      // Set X = X + R
      magma_zaxpy(mm, magma_one, d_r, 1, d_a, 1, queue);
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
  magma_free(d_a_orig);
  magma_free(d_ax);
  magma_free(d_r);
  magma_free(d_id_diag);
  return err;
}

#endif
