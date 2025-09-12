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
 */

#ifndef INCLUDE_ALGEBRAIC_H_
#define INCLUDE_ALGEBRAIC_H_

using namespace std;

/*! \brief Perform in-place matrix inversion.
 *
 * \param mat: `complex double **` The matrix to be inverted (must be a square matrix).
 * \param size: `np_int` The size of the matrix (i.e. the number of its rows or columns).
 * \param ier: `int &` Reference to an integer variable for returning a result flag.
 * \param maxrefiters: `int &` Reference to the maximum number of refinement iterations.
 * \param accuracygoal: `double &` Reference to the requested accuracy level.
 * \param refinemode: `int` Flag for refinement mode selection.
 * \param output_path: `const string &` Path where the output needs to be placed.
 * \param jxi488: `int` Index of the current wavelength calculation.
 * \param max_size: `np_int` The maximum expected size (required by some call-backs, optional, defaults to 0).
 * \param target_device: `int` ID of target GPU, if available (defaults to 0).
 */
void invert_matrix(dcomplex **mat, np_int size, int &ier, int &maxrefiters, double &accuracygoal, int refinemode, const string& output_path, int jxi488, np_int max_size=0, int target_device=0);

#endif
