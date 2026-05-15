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

/**
 * \file utils.h
 *
 * \brief Definition of auxiliary code utilities.
 */

#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

/**
 * \brief Obtain an estimate of the RAM overhead factor for optimization.
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

/**
 * \brief Write a double complex matrix to a text file.
 *
 * \param af: `VirtualAsciiFile *` Pointer to an existing VirtualAsciiFile. [IN]
 * \param mat: `dcomplex **` Pointer to the matrix. [IN]
 * \param rows: `int` Number of rows in the matrix. [IN]
 * \param columns: `int` Number of columns in the matrix. [IN]
 * \param format: `const string&` Format of the line (default is \" %5d %5d (%17.8lE,%17.8lE)\n\"). [IN]
 * \param first_index: `int` Index of the first element (default is 1, i.e. base 1 FORTRAN array notation). [IN]
 * \return result: `int` An exit code (0 if successful).
 */
int write_dcomplex_matrix(
  VirtualAsciiFile *af, dcomplex **mat, int rows, int columns,
  const std::string& format=" %5d %5d (%17.8lE,%17.8lE)\n", int first_index=1
);

/**
 * \brief Draw a PPM representation of a matrix.
 *
 * \param[in] A: `const dcomplex *` Vector representation of the matrix.
 * \param[in] m: `int64_t` Number of rows in the matrix.
 * \param[in] n: `int64_t` Number of columns in the matrix.
 * \param[in] file_name: `const string&` Reference to a string with the name of the output file.
 * \param[in] mode: `const string&` Reference to a string for the value to be represented ("MAG" to
 * draw magnitudes, "RE" to draw real parts, "IM" to draw imaginary parts. Optional, default is
 * "MAG").
 * \param[in] row_bin: `int` Number of rows to bin together (optional, default is 1).
 * \param[in] col_bin: `int` Number of columns to bin together (optional, default is 1).
 * \return result: `int` An exit code (0 if successful).
 */
int write_matrix_as_ppm(
  dcomplex *A, int64_t m, int64_t n, const std::string& file_name, const std::string& mode="MAG",
  int row_bin=1, int col_bin=1
);
#endif
