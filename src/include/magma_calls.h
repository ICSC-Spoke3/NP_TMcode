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

/*! \file magma_calls.h
 *
 * \brief C++ interface to MAGMA calls.
 *
 */

#include <string>

#ifndef INCLUDE_MAGMA_CALLS_H_
#define INCLUDE_MAGMA_CALLS_H_

/*! \brief Invert a complex matrix with double precision elements.
 *
 * call magma_zinvert1() to do the actual inversion
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 */
void magma_zinvert(dcomplex **mat, np_int n, int &jer, int device_id=0);

/*! \brief Auxiliary function to Invert a complex matrix with double precision elements.
 *
 * Use MAGMA to perform an in-place matrix inversion for a complex
 * matrix with double precision elements.
 *
 * \param inva: reference to the pointer to the first element of the matrix to be inverted on entry, inverted matrix on exit.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 * \param device_id: `int` ID of the device for matrix inversion offloading.
 */
void magma_zinvert1(dcomplex * &inva, np_int n, int &jer, int device_id);

/*! \brief Invert a complex matrix with double precision elements, applying iterative refinement of the solution
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
void magma_zinvert_and_refine(dcomplex **mat, np_int n, int &jer, int &maxrefiters, double &accuracygoal, int refinemode, int device_id, const string& output_path, int jxi488);

/*! \brief Apply iterative refinement of the solution of a matrix inversion.
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
void magma_refine(dcomplex *aorig, dcomplex *inva, np_int n, int &jer, int &maxrefiters, double &accuracygoal, int refinemode, int device_id, const string& output_path, int jxi488);

#endif
