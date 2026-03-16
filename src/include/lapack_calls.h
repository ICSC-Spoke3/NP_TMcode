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

#ifndef LAPACK_SUCCESS
#define LAPACK_SUCCESS 0
#endif

#ifdef USE_MKL
typedef MKL_Complex16 lcomplex;
#else
typedef dcomplex lcomplex;
#endif // USE_MKL

/**
 * \brief Invert a complex matrix with double precision elements.
 *
 * Use LAPACKE64 to perform an in-place matrix inversion for a complex
 * matrix with double precision elements.
 *
 * \param mat: Matrix of complex. The matrix to be inverted.
 * \param n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param jer: `int &` Reference to an integer return flag.
 * \param rs: `const RuntimeSettings &` Runtime settings instance.
 */
void zinvert(dcomplex **mat, np_int n, int &jer, const RuntimeSettings& rs=RuntimeSettings());

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
 * \param a_orig: `lcomplex *` Pointer to the first element of the non-inverted matrix. [IN]
 * \param m: `const lapack_int` Number of rows / columns in a. [IN]
 * \param a: `lcomplex *` Pointer to the inverted matrix. [IN/OUT]
 * \return err: `lapack_int` An error code (LAPACK_SUCCESS, if everything was fine).
 */
lapack_int lapack_newton(
  const RuntimeSettings& rs, lcomplex* a_orig, const lapack_int m, lcomplex* a
);

#endif
