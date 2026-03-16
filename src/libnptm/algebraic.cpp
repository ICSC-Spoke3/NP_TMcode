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

/*! \file algebraic.cpp
 *
 * \brief Implementation of algebraic functions with different call-backs.
 */

#include <string>
using namespace std;

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifdef USE_LAPACK
#ifndef INCLUDE_LAPACK_CALLS_H_
#include "../include/lapack_calls.h"
#endif
#endif

#ifdef USE_MAGMA
#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif
#endif

#ifdef USE_CUBLAS
#ifndef INCLUDE_CUBLAS_CALLS_H_
#include "../include/cublas_calls.h"
#endif
#endif

#ifndef INCLUDE_ALGEBRAIC_H_
#include "../include/algebraic.h"
#endif

// >>> FALL-BACK FUNCTIONS DECLARATION <<< //
extern void lucin(dcomplex **am, const np_int nddmst, np_int n, int &ier);
// >>>   END OF FALL-BACK FUNCTIONS    <<< //

using namespace std;

void invert_matrix(
  dcomplex **mat, np_int size, int &ier, const string& output_path, int jxi488, np_int max_size,
  int target_device, const RuntimeSettings& rs
) {
  ier = 0;
#ifdef USE_MAGMA
  magma_zinvert(mat, size, ier, target_device, rs);
#elif defined USE_CUBLAS
  cublas_zinvert(mat, size, target_device, rs);
#elif defined USE_LAPACK
  zinvert(mat, size, ier, rs);
#else
  lucin(mat, size, size, ier);
#endif
}
