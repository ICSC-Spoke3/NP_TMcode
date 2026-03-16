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

/*! \file cublas_calls.h
 *
 * \brief C++ interface to CUBLAS calls.
 *
 */

#ifndef INCLUDE_CUBLAS_CALLS_H_
#define INCLUDE_CUBLAS_CALLS_H_

/*! \brief Invert a complex matrix with double precision elements.
 *
 * Use CUBLAS to perform an in-place matrix inversion for a complex
 * matrix with double precision elements.
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 */
void cublas_zinvert(dcomplex **mat, np_int n, int device_id);

/*! \brief Invert a complex matrix with double precision elements, applying iterative refinement of the solution
 *
 * Use CUBLAS to perform matrix inversion for a complex
 * matrix with double precision elements.
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param maxrefiters: `int` Maximum number of refinement iterations to apply.
 * \param accuracygoal: `double &` Accuracy to achieve in iterative refinement, defined as the module of the maximum difference between the identity matrix and the matrix product of the (approximate) inverse times the original matrix. On return, it contains the actually achieved accuracy.
 * \param refinemode: `int` Flag to choose the refinement mode.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 */
void cublas_zinvert_and_refine(dcomplex **mat, np_int n, int &maxrefiters, double &accuracygoal, int refinemode, int device_id);

#endif
