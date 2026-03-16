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
 * \file magma_calls.h
 *
 * \brief C++ interface to MAGMA calls.
 *
 */

#include <string>

#ifndef INCLUDE_MAGMA_CALLS_H_
#define INCLUDE_MAGMA_CALLS_H_

/**
 * \brief Invert a complex matrix with double precision elements.
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 * \param rs: `const RuntimeSettings &` Runtime settings instance.
 */
void magma_zinvert(
  dcomplex **mat, np_int n, int &jer, int device_id=0,
  const RuntimeSettings& rs=RuntimeSettings()
);

/**
 * \brief Perform Newton-Schulz iterative refinement of matrix inversion.
 *
 * In this function the residual of the inversion of a matrix A is evaluated as:
 *
 * R = A^-1 A - I
 *
 * and the convergence of refinement is estimated through the largest element
 * modulus left in R.
 *
 * \param rs: `const RuntimeSettings &` Runtime settings instance. [IN]
 * \param a: `magmaDoubleComplex *` Pointer to the first element of the non-inverted matrix on host. [IN]
 * \param m: `const magma_int_t` Number of rows / columns in a. [IN]
 * \param d_a: `magmaDoubleComplex *` Pointer to the matrix on the GPU. [IN/OUT]
 * \param queue: `magma_queue_t` GPU communication queue. [IN]
 * \return err: `magma_int_t` An error code (MAGMA_SUCCESS, if everything was fine).
 */
magma_int_t magma_newton(
  const RuntimeSettings& rs, magmaDoubleComplex* a, const magma_int_t m,
  magmaDoubleComplex* d_a, magma_queue_t queue
);

/**
 * \brief Perform norm-based Newton-Schulz iterative refinement of matrix inversion.
 *
 * In this function the residual of the inversion of a matrix A is evaluated as:
 *
 * R = A A^-1 A - A
 *
 * and the convergence of refinement is estimated through the norm of R.
 *
 * \param rs: `const RuntimeSettings &` Runtime settings instance. [IN]
 * \param a: `magmaDoubleComplex *` Pointer to the first element of the non-inverted matrix on host. [IN]
 * \param m: `const magma_int_t` Number of rows / columns in a. [IN]
 * \param d_a: `magmaDoubleComplex *` Pointer to the matrix on the GPU. [IN/OUT]
 * \param queue: `magma_queue_t` GPU communication queue. [IN]
 */
magma_int_t magma_newton_norm(
  const RuntimeSettings& rs, magmaDoubleComplex* a, const magma_int_t m,
  magmaDoubleComplex* d_a, magma_queue_t queue
);

/* \brief Invert a complex matrix with double precision elements, applying iterative refinement of the solution
 *
 * call magma_zinvert1() to perform the first matrix inversion, then magma_refine() to do the refinement (only if maxrefiters is >0)
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 * \param maxrefiters: `int` Maximum number of refinement iterations to apply.
 * \param accuracygoal: `double &` Accuracy to achieve in iterative refinement, defined as the module of the maximum difference between the identity matrix and the matrix product of the (approximate) inverse times the original matrix. On return, it contains the actually achieved accuracy.
 * \param refinemode: `int` Flag to control the refinement mode.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 * \param output_path: `const string &` Path where the output needs to be placed.
 * \param jxi488: `int` Index of the current wavelength calculation.
 */
//void magma_zinvert_and_refine(dcomplex **mat, np_int n, int &jer, int &maxrefiters, double &accuracygoal, int refinemode, int device_id, const string& output_path, int jxi488);

/* \brief Apply iterative refinement of the solution of a matrix inversion.
 *
 * Iteratively compute and apply a correction to the inverse `inva` of the complex
 * matrix `aorig`, for a maximum number of `maxiters` times, or until achieving a
 * maximum residual better than `accuracygoal`.
 *
 * \param aorig: pointer to the first element of the matrix of complex to be inverted.
 * \param inva: pointer to the first element of inverse.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrices.
 * \param jer: `int &` Reference to an integer return flag.
 * \param maxrefiters: `int` Maximum number of refinement iterations to apply.
 * \param accuracygoal: `double` Accuracy to achieve in iterative refinement, defined as the module of the maximum difference between the identity matrix and the matrix product of the (approximate) inverse times the original matrix. On return, it contains the actually achieved accuracy.
 * \param refinemode: `int` Flag for refinement mode selection.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 * \param output_path: `const string &` Path where the output needs to be placed.
 * \param jxi488: `int` Index of the current wavelength calculation.
 */
// void magma_refine(dcomplex *aorig, dcomplex *inva, np_int n, int &jer, int &maxrefiters, double &accuracygoal, int refinemode, int device_id, const std::string& output_path, int jxi488);

#endif
