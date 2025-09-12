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

using namespace std;

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifdef USE_MAGMA

#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif

// hand define some preprocessor flags
//#define USE_ZGESV_GPU 1
//#define USE_ZGESV_RBT 1
//#define DEBUG_REFINE 1

#ifdef DEBUG_REFINE
#include <hdf5.h>
#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif
#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif
#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif
#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif
using namespace std;
#endif

#include <limits>

void magma_zinvert(dcomplex **mat, np_int n, int &jer, int device_id) {
  dcomplex *inva = mat[0];
  magma_zinvert1(inva, n, jer, device_id);
}

void magma_zinvert1(dcomplex * &inva, np_int n, int &jer, int device_id) {
  // magma_int_t result = magma_init();
  magma_int_t err = MAGMA_SUCCESS;
  magma_queue_t queue = NULL;
  magma_device_t dev = (magma_device_t)device_id;
  magma_queue_create(dev, &queue);
  magma_int_t info;
  magma_int_t m = (magma_int_t) n; // changed rows; a - mxm matrix
  magma_int_t mm = m * m; // size of a
  magmaDoubleComplex * a = (magmaDoubleComplex *) inva; // pointer to first element on host
#ifndef USE_ZGESV_RBT
  magma_int_t *piv;  // array of pivot indices
  piv = new magma_int_t[m]; // host mem.
#endif
  magmaDoubleComplex *d_id; // pointer to the diagonal of identity matrix
  magmaDoubleComplex *d_a; // pointer to first element on device
  magmaDoubleComplex *h_id;
  magmaDoubleComplex *dwork; // workspace
  magma_int_t ldwork; // size of dwork
  magmaDoubleComplex magma_zero;
  magma_zero.x = magma_zero.y = 0;
  magmaDoubleComplex magma_one;
  magma_one.x = 1;
  magma_one.y = 0;
  magmaDoubleComplex magma_mone;
  magma_mone.x = -1;
  magma_mone.y = 0;
#if defined(USE_ZGESV_GPU) || defined(USE_ZGESV_RBT)
  // apparently magma does not like me to set the whole identity matrix on the device from scalars with incx=0. So do set the whole damn matrix on the host, copy it on the device, destroy it on the host
  h_id = new magmaDoubleComplex[mm]();
  for (np_int i=0; i<mm; i++) {
    h_id[i] = magma_zero;
  }
  for (np_int i=0; i<m; i++) {
    h_id[i*(m+1)] = magma_one;
  }
#else
  ldwork = m * magma_get_zgetri_nb(m); // optimal block size
  err = magma_zmalloc(&dwork, ldwork); // dev. mem. for ldwork
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating dwork\n");
    exit(1);
  }
#endif
#ifdef USE_ZGESV_GPU
  err = magma_zmalloc(&d_id, mm);
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_id\n");
    exit(1);
  }
  magma_zsetmatrix(m, m, h_id, m, d_id , m, queue); // copy identity matrix to device
  delete[] h_id;
#endif
#ifndef USE_ZGESV_RBT
  // allocate matrix on device
  err = magma_zmalloc(&d_a, mm); // device memory for a
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_a\n");
    exit(1);
  }
  magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
#endif
#ifdef USE_ZGESV_GPU
  magma_zgesv_gpu(m, m, d_a, m, piv, d_id, m, &info);
  delete[] piv; // free host memory
  magma_free(d_a);
  magma_zgetmatrix(m, m, d_id , m, a, m, queue); // copy d_id -> a
  magma_free(d_id);
#elif defined USE_ZGESV_RBT
  magma_zgesv_rbt(MagmaTrue, m, m, a, m, h_id, m, &info);
  memcpy(a, h_id, mm*sizeof(magmaDoubleComplex));
  delete[] h_id;
#else
  magma_zgetrf_gpu(m, m, d_a, m, piv, &info);
  magma_zgetri_gpu(m, d_a, m, piv, dwork, ldwork, &info);
  magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy d_a -> a
  delete[] piv; // free host memory
  magma_free(d_a);
  magma_free(dwork);
#endif
  magma_queue_destroy(queue); // destroy queue
  // result = magma_finalize();
  jer = (int)err;
}

void magma_refine(dcomplex *aorig, dcomplex *inva, np_int n, int &jer, int &maxiters, double &accuracygoal, int refinemode, int device_id, const string& output_path, int jxi488) {
  magma_int_t err = MAGMA_SUCCESS;
  magma_queue_t queue = NULL;
  magma_device_t dev = (magma_device_t)device_id;
  magma_queue_create(dev, &queue);
  magma_int_t info;
  magma_int_t m = (magma_int_t)n; // changed rows; a - mxm matrix
  magma_int_t mm = m * m; // size of a
  magmaDoubleComplex *a = (magmaDoubleComplex *) inva;
  magmaDoubleComplex *a_orig = (magmaDoubleComplex *) aorig;
  magmaDoubleComplex *d_a; // pointer to first element on device
  magmaDoubleComplex *d_a_orig; // pointer to original array on device
  magmaDoubleComplex *d_a_residual; // pointer to residual array on device
  magmaDoubleComplex *d_a_refine; // pointer to residual array on device
  magmaDoubleComplex *d_id0; // pointer to the diagonal of identity matrix
  magmaDoubleComplex magma_zero;
  magma_zero.x = magma_zero.y = 0;
  magmaDoubleComplex magma_one;
  magma_one.x = 1;
  magma_one.y = 0;
  magmaDoubleComplex magma_mone;
  magma_mone.x = -1;
  magma_mone.y = 0;
#ifdef DEBUG_REFINE
  char virtual_line[256];
  dcomplex **ambuf = new dcomplex * [n]() ;
  ambuf[0] = new dcomplex[n*n]();
  for (np_int i=1; i<n; i++) {
    ambuf[i] = ambuf[0] + i*n ;
  }
  magmaDoubleComplex *abuf = (magmaDoubleComplex *) ambuf[0];
  memcpy(ambuf[0], a_orig, n*n*sizeof(dcomplex));
  VirtualAsciiFile *outam = new VirtualAsciiFile();
  string outam_name = output_path + "/c_AM4_JXI" + to_string(jxi488) + ".txt";
  sprintf(virtual_line, " AM matrix before inversion\n");
  outam->append_line(virtual_line);
#ifdef USE_ILP64
  sprintf(virtual_line, " %ld\n", n);
#else
  sprintf(virtual_line, " %d\n", n);
#endif
  outam->append_line(virtual_line);  
  sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
  outam->append_line(virtual_line);
  write_dcomplex_matrix(outam, ambuf, n, n);
  outam->write_to_disk(outam_name);
  delete outam;
  memcpy(ambuf[0], a, n*n*sizeof(dcomplex));
  outam = new VirtualAsciiFile();
  outam_name = output_path + "/c_AM5_JXI" + to_string(jxi488) + "_iter0.txt";
  sprintf(virtual_line, " AM inverse before refinement\n");
  outam->append_line(virtual_line);
#ifdef USE_ILP64
  sprintf(virtual_line, " %ld\n", n);
#else
  sprintf(virtual_line, " %d\n", n);
#endif
  outam->append_line(virtual_line);  
  sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
  outam->append_line(virtual_line);
  write_dcomplex_matrix(outam, ambuf, n, n);
  outam->write_to_disk(outam_name);
  delete outam;
#endif
  err = magma_zmalloc(&d_id0, 1);
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_id0\n");
    exit(1);
  }
  magma_zsetvector(1, &magma_one, 1, d_id0, 1, queue);
  err = magma_zmalloc(&d_a_orig, mm); // device memory for the original a
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_a_orig\n");
    exit(1);
  }
  magma_zsetmatrix(m, m, a_orig, m, d_a_orig , m, queue);
  err = magma_zmalloc(&d_a, mm); // device memory for the inverse matrix
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_a\n");
    exit(1);
  }
  magma_zsetmatrix(m, m, a, m, d_a , m, queue);
  err = magma_zmalloc(&d_a_residual, mm); // device memory for iterative correction of inverse of a
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_a_residual\n");
    exit(1);
  }
  err = magma_zmalloc(&d_a_refine, mm); // device memory for iterative correction of inverse of a
  if (err != MAGMA_SUCCESS) {
    printf("Error allocating d_a_refine\n");
    exit(1);
  }
  bool iteraterefine = true;
  // multiply minus the original matrix times the inverse matrix
  // NOTE: factors in zgemm are swapped because zgemm is designed for column-major
  // Fortran-style arrays, whereas our arrays are C-style row-major.
  magma_zgemm(MagmaNoTrans, MagmaNoTrans, m, m, m,  magma_mone, d_a, m, d_a_orig, m, magma_zero, d_a_residual, m, queue);
  // add the identity to the product
  magma_zaxpy(m, magma_one, d_id0, 0, d_a_residual, m+1, queue);
#ifdef DEBUG_REFINE
  magma_zgetmatrix(m, m, d_a_residual, m, abuf, m, queue); // copy residual to buffer
  outam = new VirtualAsciiFile();
  outam_name = output_path + "/c_AM6_JXI" + to_string(jxi488) + "_iter0.txt";
  sprintf(virtual_line, " AM residual before refinement\n");
  outam->append_line(virtual_line);
#ifdef USE_ILP64
  sprintf(virtual_line, " %ld\n", n);
#else
  sprintf(virtual_line, " %d\n", n);
#endif
  outam->append_line(virtual_line);  
  sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
  outam->append_line(virtual_line);
  write_dcomplex_matrix(outam, ambuf, n, n);
  outam->write_to_disk(outam_name);
  delete outam;
#endif
  double oldmax=0;
  if (refinemode >0) {
    // find the maximum absolute value of the residual
    magma_int_t maxindex = magma_izamax(mm, d_a_residual, 1, queue);
    magmaDoubleComplex magmamax;
    // transfer the maximum value to the host
    magma_zgetvector(1, d_a_residual+maxindex, 1, &magmamax, 1, queue);
    // take the module
    oldmax = cabs(magmamax.x + I*magmamax.y);
    printf("Initial max residue = %g\n", oldmax);
    if (oldmax < accuracygoal) iteraterefine = false;
  }
  // begin correction loop (should iterate maxiters times)
  int iter;
  for (iter=0; (iter<maxiters) && iteraterefine; iter++) {
    // multiply the inverse times the residual, add to the initial inverse
    magma_zgemm(MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_a_residual, m, d_a, m, magma_zero, d_a_refine, m, queue);
#ifdef DEBUG_REFINE
    magma_zgetmatrix(m, m, d_a_refine, m, abuf, m, queue); // copy refinement to buffer
    outam = new VirtualAsciiFile();
    outam_name = output_path + "/c_AM7_JXI" + to_string(jxi488) + "_iter" + to_string(iter+1) + ".txt";
    sprintf(virtual_line, " AM correction refinement %d\n", iter+1);
    outam->append_line(virtual_line);
#ifdef USE_ILP64
    sprintf(virtual_line, " %ld\n", n);
#else
    sprintf(virtual_line, " %d\n", n);
#endif
    outam->append_line(virtual_line);  
    sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
    outam->append_line(virtual_line);
    write_dcomplex_matrix(outam, ambuf, n, n);
    outam->write_to_disk(outam_name);
    delete outam;
#endif
    // add to the initial inverse
    magma_zaxpy (mm, magma_one, d_a_refine, 1, d_a, 1, queue);
#ifdef DEBUG_REFINE
    magma_zgetmatrix(m, m, d_a, m, abuf, m, queue); // copy new inverse to buffer
    outam = new VirtualAsciiFile();
    outam_name = output_path + "/c_AM5_JXI" + to_string(jxi488) + "_iter" + to_string(iter+1) + ".txt";
    sprintf(virtual_line, " AM inverse after refinement %d\n", iter+1);
    outam->append_line(virtual_line);
#ifdef USE_ILP64
    sprintf(virtual_line, " %ld\n", n);
#else
    sprintf(virtual_line, " %d\n", n);
#endif
    outam->append_line(virtual_line);  
    sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
    outam->append_line(virtual_line);
    write_dcomplex_matrix(outam, ambuf, n, n);
    outam->write_to_disk(outam_name);
    delete outam;
#endif
    // multiply minus the original matrix times the new inverse matrix
    magma_zgemm(MagmaNoTrans, MagmaNoTrans, m, m, m, magma_mone, d_a, m, d_a_orig, m, magma_zero, d_a_residual, m, queue);
    // add the identity to the product
    magma_zaxpy (m, magma_one, d_id0, 0, d_a_residual, m+1, queue);
#ifdef DEBUG_REFINE
    magma_zgetmatrix(m, m, d_a_residual, m, abuf, m, queue); // copy residual to buffer
    outam = new VirtualAsciiFile();
    outam_name = output_path + "/c_AM6_JXI" + to_string(jxi488) + "_iter" + to_string(iter+1) + ".txt";
    sprintf(virtual_line, " AM residual after refinement %d\n", iter+1);
    outam->append_line(virtual_line);
#ifdef USE_ILP64
    sprintf(virtual_line, " %ld\n", n);
#else
    sprintf(virtual_line, " %d\n", n);
#endif
    outam->append_line(virtual_line);  
    sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
    outam->append_line(virtual_line);
    write_dcomplex_matrix(outam, ambuf, n, n);
    outam->write_to_disk(outam_name);
    delete outam;
#endif
    if ((refinemode==2) || ((refinemode==1) && (iter == (maxiters-1)))) {
      // find the maximum absolute value of the residual
      magma_int_t maxindex = magma_izamax(mm, d_a_residual, 1, queue);
      // transfer the maximum value to the host
      magmaDoubleComplex magmamax;
      magma_zgetvector(1, d_a_residual+maxindex, 1, &magmamax, 1, queue);
      // take the module
      double newmax = cabs(magmamax.x + I*magmamax.y);
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
  // free temporary device arrays 
  magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy final refined d_a -> a
#ifdef DEBUG_REFINE
  delete[] ambuf[0];
  delete[] ambuf;
#endif
  magma_free(d_a); // free device memory
  magma_free(d_id0);
  magma_free(d_a_orig); // free device memory
  magma_free(d_a_residual); // free device memory
  magma_free(d_a_refine); // free device memory
  magma_queue_destroy(queue); // destroy queue
  jer = (int)err;
}

void magma_zinvert_and_refine(dcomplex **mat, np_int n, int &jer, int &maxiters, double &accuracygoal, int refinemode, int device_id, const string& output_path, int jxi488) {
  if (maxiters>0) {
    dcomplex *inva = new dcomplex[n*n]();
    dcomplex *aorig = mat[0];
    memcpy(inva, aorig, n*n*sizeof(dcomplex));
    magma_zinvert1(inva, n, jer, device_id);
    magma_refine(aorig, inva, n, jer, maxiters, accuracygoal, refinemode, device_id, output_path, jxi488);
    memcpy(aorig, inva, n*n*sizeof(dcomplex));
    delete[] inva;
  }
  else {
    dcomplex *inva = mat[0];
    magma_zinvert1(inva, n, jer, device_id);
  }
}
    
//   // magma_int_t result = magma_init();
//   magma_int_t err = MAGMA_SUCCESS;
//   magma_queue_t queue = NULL;
//   magma_device_t dev = (magma_device_t)device_id;
//   magma_queue_create(dev, &queue);
//   magma_int_t info;
//   magma_int_t m = (magma_int_t)n; // changed rows; a - mxm matrix
//   magma_int_t mm = m * m; // size of a
// #ifndef USE_ZGESV_RBT
//   magma_int_t *piv;  // array of pivot indices
//   piv = new magma_int_t[m]; // host mem.
// #endif
//   magmaDoubleComplex *a = (magmaDoubleComplex *) mat[0]; // pointer to first element on host
// #ifdef USE_ZGESV_RBT
//   magmaDoubleComplex * a_orig;
// #endif
//   magmaDoubleComplex *d_a; // pointer to first element on device
//   magmaDoubleComplex *d_a_orig; // pointer to original array on device
//   magmaDoubleComplex *d_a_residual; // pointer to residual array on device
//   magmaDoubleComplex *d_a_refine; // pointer to residual array on device
//   magmaDoubleComplex *d_id0; // pointer to the diagonal of identity matrix
// #ifdef USE_ZGESV_GPU
//   magmaDoubleComplex *d_id; // pointer to the diagonal of identity matrix
// #endif
//   magmaDoubleComplex native_zero;
//   native_zero.x = native_zero.y = 0;
//   magmaDoubleComplex native_one;
//   native_one.x = 1;
//   native_one.y = 0;
//   if (maxiters>0) {
//     err = magma_zmalloc(&d_id0, 1);
//     if (err != MAGMA_SUCCESS) {
//       printf("Error allocating d_id0\n");
//       exit(1);
//     }
//     magma_zsetvector(1, &native_one, 1, d_id0, 1, queue);
//   }
// #if defined(USE_ZGESV_GPU) || defined(USE_ZGESV_RBT)
//   // apparently magma does not like me to set the whole identity matrix on the device from scalars with incx=0. So do set the whole damn matrix on the host, copy it on the device, destroy it on the host
//   magmaDoubleComplex *h_id = new magmaDoubleComplex[mm]();
//   for (np_int i=0; i<mm; i++) {
//     h_id[i] = native_zero;
//   }
//   for (np_int i=0; i<m; i++) {
//     h_id[i*(m+1)] = native_one;
//   }
// #else
//   magmaDoubleComplex *dwork; // workspace
//   magma_int_t ldwork; // size of dwork
//   ldwork = m * magma_get_zgetri_nb(m); // optimal block size
//   err = magma_zmalloc(&dwork, ldwork); // dev. mem. for ldwork
//   if (err != MAGMA_SUCCESS) {
//     printf("Error allocating dwork\n");
//     exit(1);
//   }
// #endif
// #ifdef USE_ZGESV_GPU
//   err = magma_zmalloc(&d_id, mm);
//   if (err != MAGMA_SUCCESS) {
//     printf("Error allocating d_id\n");
//     exit(1);
//   }
//   magma_zsetmatrix(m, m, h_id, m, d_id , m, queue); // copy identity matrix to device
//   delete[] h_id;
// #endif
// #ifdef USE_ZGESV_RBT
//   if (maxiters>0) {
//     a_orig =  new magmaDoubleComplex[mm]();
//     memcpy(a_orig, a, mm*sizeof(magmaDoubleComplex));
//     // for (np_int i=0; i<mm; i++) {
//     //   a_orig[i] = a[i];
//     // }
//   }
//   magma_zgesv_rbt(MagmaTrue, m, m, a, m, h_id, m, &info);
//   err = magma_zmalloc(&d_a, mm); // device memory for a, copy the inverse to it
//   if (err != MAGMA_SUCCESS) {
//     printf("Error allocating d_a\n");
//     exit(1);
//   }
//   magma_zsetmatrix(m, m, h_id, m, d_a , m, queue);
//   delete[] h_id;
//   if (maxiters>0) {
//     err = magma_zmalloc(&d_a_orig, mm); // device memory for copy of a
//     if (err != MAGMA_SUCCESS) {
//       printf("Error allocating d_a_orig\n");
//       exit(1);
//     }
//     magma_zsetmatrix(m, m, a_orig, m, d_a_orig , m, queue);
//     delete[] a_orig;
//   }
// #elif defined USE_ZGESV_GPU
//   err = magma_zmalloc(&d_a, mm); // device memory for a
//   if (err != MAGMA_SUCCESS) {
//     printf("Error allocating d_a\n");
//     exit(1);
//   }
//   magma_zsetmatrix(m, m, a, m, d_a , m, queue);
//   magma_zgesv_gpu(m, m, d_a, m, piv, d_id, m, &info);
//   if (maxiters>0) {
//     d_a_orig = d_a;
//     // copy again the original matrix on the device, in d_a_orig
//     magma_zsetmatrix(m, m, a, m, d_a_orig , m, queue);
//   }
//   else {
//     magma_free(d_a);
//   }
//   d_a = d_id;
//   delete[] piv; // free host memory
//   // magma_free(d_a);
// #else
//   err = magma_zmalloc(&d_a, mm); // device memory for a
//   if (err != MAGMA_SUCCESS) {
//     printf("Error allocating d_a\n");
//     exit(1);
//   }
//   magma_zsetmatrix(m, m, a, m, d_a , m, queue);
//   magma_zgetrf_gpu(m, m, d_a, m, piv, &info);
//   magma_zgetri_gpu(m, d_a, m, piv, dwork, ldwork, &info);
//   delete[] piv; // free host memory
//   magma_free(dwork);
//   if (maxiters>0) {
//     err = magma_zmalloc(&d_a_orig, mm); // device memory for copy of a
//     if (err != MAGMA_SUCCESS) {
//       printf("Error allocating d_a_orig\n");
//       exit(1);
//     }
//     magma_zsetmatrix(m, m, a, m, d_a_orig , m, queue);
//   }
// #endif

//   if (maxiters>0) {
//     // allocate memory for the temporary matrix products
//     err = magma_zmalloc(&d_a_residual, mm); // device memory for iterative correction of inverse of a
//     if (err != MAGMA_SUCCESS) {
//       printf("Error allocating d_a_residual\n");
//       exit(1);
//     }
//     err = magma_zmalloc(&d_a_refine, mm); // device memory for iterative correction of inverse of a
//     if (err != MAGMA_SUCCESS) {
//       printf("Error allocating d_a_refine\n");
//       exit(1);
//     }
//   }
//   bool iteraterefine = true;
//   if (maxiters>0) {
//     magmaDoubleComplex magma_mone;
//     magma_mone.x = -1;
//     magma_mone.y = 0;
//     magmaDoubleComplex magma_one;
//     magma_one.x = 1;
//     magma_one.y = 0;
//     magmaDoubleComplex magma_zero;
//     magma_zero.x = 0;
//     magma_zero.y = 0;
//     // multiply minus the original matrix times the inverse matrix
//     // NOTE: factors in zgemm are swapped because zgemm is designed for column-major
//     // Fortran-style arrays, whereas our arrays are C-style row-major.
//     magma_zgemm(MagmaNoTrans, MagmaNoTrans, m, m, m,  magma_mone, d_a, m, d_a_orig, m, magma_zero, d_a_residual, m, queue);
//     // add the identity to the product
//     magma_zaxpy(m, magma_one, d_id0, 0, d_a_residual, m+1, queue);
//     double oldmax=0;
//     if (refinemode >0) {
//       // find the maximum absolute value of the residual
//       magma_int_t maxindex = magma_izamax(mm, d_a_residual, 1, queue);
//       magmaDoubleComplex magmamax;
//       // transfer the maximum value to the host
//       magma_zgetvector(1, d_a_residual+maxindex, 1, &magmamax, 1, queue);
//       // take the module
//       oldmax = cabs(magmamax.x + I*magmamax.y);
//       printf("Initial max residue = %g\n", oldmax);
//       if (oldmax < accuracygoal) iteraterefine = false;
//     }
//     // begin correction loop (should iterate maxiters times)
//     int iter;
//     for (iter=0; (iter<maxiters) && iteraterefine; iter++) {
//       // multiply the inverse times the residual, add to the initial inverse
//       magma_zgemm(MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_a_residual, m, d_a, m, magma_zero, d_a_refine, m, queue);
//       // add to the initial inverse
//       magma_zaxpy (mm, magma_one, d_a_refine, 1, d_a, 1, queue);
//       // multiply minus the original matrix times the new inverse matrix
//       magma_zgemm(MagmaNoTrans, MagmaNoTrans, m, m, m, magma_mone, d_a, m, d_a_orig, m, magma_zero, d_a_residual, m, queue);
//       // add the identity to the product
//       magma_zaxpy (m, magma_one, d_id0, 0, d_a_residual, m+1, queue);
//       if ((refinemode==2) || ((refinemode==1) && (iter == (maxiters-1)))) {
// 	// find the maximum absolute value of the residual
// 	magma_int_t maxindex = magma_izamax(mm, d_a_residual, 1, queue);
// 	// transfer the maximum value to the host
// 	magmaDoubleComplex magmamax;
// 	magma_zgetvector(1, d_a_residual+maxindex, 1, &magmamax, 1, queue);
// 	// take the module
// 	double newmax = cabs(magmamax.x + I*magmamax.y);
// 	printf("Max residue after %d iterations = %g\n", iter+1, newmax);
// 	// if the maximum in the residual decreased from the previous iteration,
// 	// update oldmax and go on, otherwise no point further iterating refinements
// 	if ((refinemode==2) && ((newmax > oldmax)||(newmax < accuracygoal))) iteraterefine = false;
// 	oldmax = newmax;
//       }
//     }
//     // if we are being called with refinemode=2, then on exit we set maxiters to the actual number of iters we performed to achieve the required accuracy
//     if (refinemode==2) maxiters = iter;
//     accuracygoal = oldmax;
//     // end correction loop
//   }
//   // free temporary device arrays 
//   magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy final refined d_a -> a
//   magma_free(d_a); // free device memory
//   // delete[] a_unref;
//   if (maxiters>0) {
//     magma_free(d_id0);
//     magma_free(d_a_orig); // free device memory
//     magma_free(d_a_residual); // free device memory
//     magma_free(d_a_refine); // free device memory
//   }
//   magma_queue_destroy(queue); // destroy queue
//   // result = magma_finalize();
//   jer = (int)err;
// }


#endif
