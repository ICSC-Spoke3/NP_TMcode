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

/*! \file types.h
 *
 * \brief Definition of fundamental types in use.
 */

#ifndef INCLUDE_TYPES_H_
#define INCLUDE_TYPES_H_

#include <complex.h>

//! \brief Short-cut to C-style double precision complex type.
typedef __complex__ double dcomplex;

#ifdef USE_MKL
#ifdef USE_ILP64
#ifndef MKL_INT 
#define MKL_INT int64_t
#endif // MKL_INT
#else
#ifndef MKL_INT 
#define MKL_INT int32_t
#endif // MKL_INT
#endif // USE_ILP64
#endif // USE_MKL

#ifdef USE_LAPACK
#ifdef USE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif // USE_MKL
#endif // USE_LAPACK

#ifdef USE_MAGMA
#include <magma_v2.h>
#endif

#ifndef np_int
#ifdef lapack_int
//! \brief Alias of the default integer type.
#define np_int lapack_int
#else
#ifdef USE_ILP64
//! \brief Alias of the default integer type.
#define np_int int64_t
#else
//! \brief Alias of the default integer type.
#define np_int int32_t
#endif // USE_ILP64
#endif // lapack_int
#endif // np_int

//! \brief Macro to compute the conjugate of a complex number.
#define dconjg(z) ( (__real__ (z) - I * (__imag__ (z))) )

/*! \brief Get the imaginary part of double precision complex number.
 *
 * \param z: `const dcomplex &` Reference to the complex number from
 * which to extract the imaginary part.
 * \return Im(z): `double` Imaginary part of z.
 */
double inline imag(const dcomplex &z) { return __imag__ z; }

/*! \brief Get the imaginary part of double precision complex number.
 *
 * \param z: `const dcomplex &` Reference to the complex number from
 * which to extract the real part.
 * \return Re(z): `double` Real part of z.
 */
double inline real(const dcomplex &z) { return __real__ z; }
#endif
