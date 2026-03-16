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

/*! \file lapack_calls.h
 *
 * \brief C++ interface to LAPACK calls.
 *
 */

#ifndef INCLUDE_LAPACK_CALLS_H_
#define INCLUDE_LAPACK_CALLS_H_

/*! \brief Invert a complex matrix with double precision elements.
 *
 * Use LAPACKE64 to perform an in-place matrix inversion for a complex
 * matrix with double precision elements.
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 */
void zinvert(dcomplex **mat, np_int n, int &jer);

/*! \brief Invert a complex matrix with double precision elements, using iterative refinement to improve accuracy.
 *
 * Use LAPACKE64 to perform an in-place matrix inversion for a complex
 * matrix with double precision elements.
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 * \param maxrefiters: `int` Maximum number of refinement iterations.
 * \param accuracygoal: `double &` Accuracy to achieve in iterative refinement, defined as the module of the maximum difference between the identity matrix and the matrix product of the (approximate) inverse times the original matrix. On return, it contains the actually achieved accuracy.
 * \param refinemode: `int` Flag to control the refinement mode.
 */
void zinvert_and_refine(dcomplex **mat, np_int n, int &jer, int &maxrefiters, double &accuracygoal, int refinemode);

#endif
