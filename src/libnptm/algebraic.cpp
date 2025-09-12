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

#ifdef USE_LAPACK
// define by hand for a first test
//#define USE_REFINEMENT 1
#ifndef INCLUDE_LAPACK_CALLS_H_
#include "../include/lapack_calls.h"
#endif
#endif

#ifdef USE_MAGMA
// define by hand for a first test
//#define USE_REFINEMENT 1
#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif
#endif

// define by hand for a first test
//#define USE_CUBLAS 1
#ifdef USE_CUBLAS
// define by hand for a first test
//#define USE_REFINEMENT 1
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

void invert_matrix(dcomplex **mat, np_int size, int &ier, int &maxrefiters, double &accuracygoal, int refinemode, const string& output_path, int jxi488, np_int max_size, int target_device) {
  ier = 0;
#ifdef USE_MAGMA
#ifdef USE_REFINEMENT
  // try using the iterative refinement to obtain a more accurate solution
  // we pass to magma_zinvert_and_refine() the accuracygoal in, get the actual
  // accuracy back out
  magma_zinvert_and_refine(mat, size, ier, maxrefiters, accuracygoal, refinemode, target_device, output_path, jxi488);
#else
  magma_zinvert(mat, size, ier, target_device);
#endif  
#elif defined USE_CUBLAS
#ifdef USE_REFINEMENT
  cublas_zinvert_and_refine(mat, size, maxrefiters, accuracygoal, refinemode, target_device);
#else
  cublas_zinvert(mat, size, target_device);
#endif
#elif defined USE_LAPACK
#ifdef USE_REFINEMENT
  zinvert_and_refine(mat, size, ier, maxrefiters, accuracygoal, refinemode);
#else
  zinvert(mat, size, ier);
#endif
#else
  lucin(mat, max_size, size, ier);
#endif
}
