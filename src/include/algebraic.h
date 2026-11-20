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

/*! \file algebraic.h
 *
 * \brief Declaration of algebraic functions with different call-backs.
 *
 * In principle, the system that runs NP_TMcode may offer various types of
 * optimized features, such as multi-core or multi-node scaling, GPU offload,
 * or external libraries. This header collects a set of functions that can
 * perform standard algebraic operations choosing the most optimized available
 * system as a call-back. If no optimization is detected, eventually the 
 * legacy serial function implementation is used as a fall-back.
 * \ingroup linalg
 */

#ifndef INCLUDE_ALGEBRAIC_H_
#define INCLUDE_ALGEBRAIC_H_

/*! \addtogroup linalg
 * @{
 */

using namespace std;

/**
 * \brief Check that the build can invert matrices with an inversion mode.
 *
 * MAGMA builds implement all the modes; cuBLAS builds implement LU only;
 * LAPACK builds implement LU and GESV; the fall-back `lucin()` always
 * inverts by LU decomposition.
 *
 * \param invert_mode: `short` The inversion mode (a `RuntimeSettings::INV_MODE_*` value).
 * \throws ParseError If the build does not implement the mode.
 */
void check_inversion_mode(short invert_mode);

/**
 * \brief Perform in-place matrix inversion.
 *
 * In builds with MAGMA or cuBLAS, the matrix is inverted on the GPU
 * `target_device`, if it is not negative. If the GPU fails, the device
 * policy of `rs` decides: `STOP` throws a `DeviceError`, `USE_CPU_FALLBACK`
 * writes a warning and inverts the matrix on the CPU (with LAPACK, if
 * available, and by LU decomposition for the modes that only MAGMA
 * implements).
 *
 * \param mat: `complex double **` The matrix to be inverted (must be a square matrix).
 * \param size: `np_int` The size of the matrix (i.e. the number of its rows or columns).
 * \param ier: `int &` Reference to an integer variable for returning a result flag.
 * \param output_path: `const string &` Path where the output needs to be placed.
 * \param jxi488: `int` Index of the current wavelength calculation.
 * \param max_size: `np_int` The maximum expected size (required by some call-backs, optional, defaults to 0).
 * \param target_device: `int` ID of target GPU, or -1 to use the CPU (defaults to 0).
 * \param rs: `const RuntimeSettings &` Runtime options instance (optional).
 * \throws DeviceError If the GPU fails and the device policy is `STOP`.
 */
void invert_matrix(dcomplex **mat, np_int size, int &ier, const string& output_path, int jxi488, np_int max_size=0, int target_device=0, const RuntimeSettings& rs=RuntimeSettings());

//! @}
#endif
