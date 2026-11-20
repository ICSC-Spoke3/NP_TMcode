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
 * \ingroup linalg
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
#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif
#endif

#ifdef USE_CUBLAS
#ifndef INCLUDE_CUBLAS_CALLS_H_
#include "../include/cublas_calls.h"
#endif
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_RESOURCES_H_
#include "../include/resources.h"
#endif

#ifndef INCLUDE_ALGEBRAIC_H_
#include "../include/algebraic.h"
#endif

// >>> FALL-BACK FUNCTIONS DECLARATION <<< //
extern void lucin(dcomplex **am, const np_int nddmst, np_int n, int &ier);
// >>>   END OF FALL-BACK FUNCTIONS    <<< //

using namespace std;

void check_inversion_mode(short invert_mode) {
  const char *names[] = {"LU", "GESV", "RBT", "SVD"};
  string name = (invert_mode >= 0 && invert_mode <= 3) ? names[invert_mode] : to_string(invert_mode);
#ifdef USE_MAGMA
  bool available = (invert_mode >= RuntimeSettings::INV_MODE_LU && invert_mode <= RuntimeSettings::INV_MODE_SVD);
  string build = "MAGMA";
#elif defined USE_CUBLAS
  bool available = (invert_mode == RuntimeSettings::INV_MODE_LU);
  string build = "cuBLAS (LU only)";
#elif defined USE_LAPACK
  bool available = (invert_mode == RuntimeSettings::INV_MODE_LU || invert_mode == RuntimeSettings::INV_MODE_GESV);
  string build = "LAPACK (LU and GESV only)";
#else
  bool available = true; // lucin() inverts by LU decomposition in all modes.
  string build = "";
#endif
  if (!available) {
    throw ParseError("inversion mode " + name + " is not available in this build, which uses " + build);
  }
}

//! \cond
// Invert a matrix on the CPU: with LAPACK, if available (by LU decomposition
// for the modes that LAPACK does not implement), with lucin() otherwise.
static void cpu_invert(dcomplex **mat, np_int size, int &ier, const RuntimeSettings& rs) {
#ifdef USE_LAPACK
  short mode = rs.invert_mode;
  if (mode != RuntimeSettings::INV_MODE_LU && mode != RuntimeSettings::INV_MODE_GESV) {
    mode = RuntimeSettings::INV_MODE_LU;
  }
  zinvert(mat, size, ier, rs, mode);
#else
  lucin(mat, size, size, ier);
#endif
}
//! \endcond

void invert_matrix(
  dcomplex **mat, np_int size, int &ier, const string& output_path, int jxi488, np_int max_size,
  int target_device, const RuntimeSettings& rs
) {
  ier = 0;
#if defined(USE_MAGMA) || defined(USE_CUBLAS)
  if (target_device >= 0) {
    try {
#ifdef USE_MAGMA
      magma_zinvert(mat, size, ier, target_device, rs);
#else
      cublas_zinvert(mat, size, ier, target_device, rs);
#endif
      return;
    } catch (const DeviceError &ex) {
      // The GPU wrappers leave the matrix unchanged when they fail.
      if (rs.device_policy != DEVICE_POLICY_USE_CPU_FALLBACK) throw;
      if (rs.logger != NULL) {
	rs.logger->log(
	  "GPU " + to_string(target_device) + " failed on scale " + to_string(jxi488) + " ("
	  + ex.what() + "): inverting the matrix on the CPU (DEVICE_POLICY=USE_CPU_FALLBACK).\n",
	  LOG_WARN
	);
      }
      ier = 0;
    }
  }
#endif
  cpu_invert(mat, size, ier, rs);
}
