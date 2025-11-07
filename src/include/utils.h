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

/*! \file utils.h
 *
 * \brief Definition of auxiliary code utilities.
 */

#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

/*! \brief Obtain an estimate of the RAM overhead factor for optimization.
 *
 * Code speed-up optimization usually comes at the cost of increased memory
 * requirements. While this should not generally be a serious issue, it can
 * become such in case of models that require large amounts of memory to be
 * handled. This function tests the code build configuration to provide an
 * zero-order estimate of the weight of these overheads to protect the host
 * system from saturating the RAM.
 *
 * \return factor: `double` The multiplicative factor to be applied to data size.
 */
double get_ram_overhead();

/*! \brief Write a double complex matrix to a text file.
 *
 * \param af: `VirtualAsciiFile *` Pointer to an existing VirtualAsciiFile.
 * \param mat: `dcomplex **` Pointer to the matrix.
 * \param rows: `int` Number of rows in the matrix.
 * \param columns: `int` Number of columns in the matrix.
 * \param format: `const string&` Format of the line (default is \" %5d %5d (%17.8lE,%17.8lE)\n\")
 * \param first_index: `int` Index of the first element (default is 1, i.e. base 1 FORTRAN array notation)
 */
int write_dcomplex_matrix(
			  VirtualAsciiFile *af, dcomplex **mat, int rows,
			  int columns, const std::string& format=" %5d %5d (%17.8lE,%17.8lE)\n",
			  int first_index=1
);

#endif
