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

/*! \file utils.cpp
 *
 * \brief Implementation of auxiliary code utilities.
 */

#include <hdf5.h>
#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif

using namespace std;

double get_ram_overhead() {
  double result = 1.0;
#ifdef USE_MAGMA
  if (result < 15.0) result = 15.0;
#endif //USE_MAGMA
#ifdef USE_LAPACK
  if (result < 2.0) result = 2.0;
#endif //USE_MAGMA
  return result;
}


int write_dcomplex_matrix(
  VirtualAsciiFile *af, dcomplex **mat, int rows, int columns,
  const std::string& format, int first_index
) {
  int result = 0;
  char virtual_line[256];
  for (int i=0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      sprintf(
	      virtual_line, format.c_str(), i + first_index, j + first_index,
	      real(mat[i][j]), imag(mat[i][j])
      );
      af->append_line(virtual_line);
    }
  }
  return result;
}
