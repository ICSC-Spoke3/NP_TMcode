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

#include <fstream>
#include <cstring>

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

int write_matrix_as_ppm(
  dcomplex *A, int64_t m, int64_t n, const std::string& file_name, const std::string& mode,
  int row_bin, int col_bin
) {
  int result = 0;
  double (*cabs)(const dcomplex& z) = [] (const dcomplex& z) -> double { return dcabs(z); };
  double (*function)(const dcomplex& z) = cabs;
  if (mode.compare("REAL") == 0) {
    function = &real;
  } else if (mode.compare("IMAG") == 0) {
    function = &imag;
  }
  int64_t mn = m * n;
  int64_t ppm_height = (row_bin > 1) ? m / row_bin : m;
  int64_t ppm_width = (col_bin > 1) ? n / col_bin : n;
  const int height_spare = m % row_bin;
  const int width_spare = n % col_bin;
  if (height_spare > 0) ++ppm_height;
  if (width_spare > 0) ++ppm_width;
  int64_t ppm_size = ppm_width * ppm_height;
  int bin_size = row_bin * col_bin;
  
  ofstream ofs(file_name, ios::binary);
  // Header PPM: P5 = binary grayscale, width, height, max value
  ofs << "P5\n" << ppm_width << " " << ppm_height << "\n255\n";

  double max_val = 0.0, min_val = 1.0e+09;
  double *logs = new double[ppm_size];

#pragma omp parallel for reduction(max: max_val) \
  reduction(min: min_val)
  for (int64_t pi = 0; pi < ppm_size; ++pi) {
    int this_bin_size = bin_size;
    int64_t vi = pi / ppm_width;
    int64_t vj = pi % ppm_width;
    int64_t first_i = vi * row_bin;
    int64_t last_i = first_i + row_bin;
    if (last_i > m) {
      last_i = m;
      this_bin_size /= row_bin;
      this_bin_size *= height_spare;
    }
    int64_t first_j = vj * col_bin;
    int64_t last_j = first_j + col_bin;
    if (last_j > n) {
      last_j = n;
      this_bin_size /= col_bin;
      this_bin_size *= width_spare;
    }
    double value = 0.0;
    int64_t asize = (last_i - first_i) * (last_j - first_j);
    for (int64_t aij = 0; aij < asize; ++aij) {
      int64_t ai = first_i + aij / (last_j - first_j);
      int64_t aj = first_j + aij % (last_j - first_j);
      value += function(A[n * ai + aj]);
    }
    value /= this_bin_size;
    logs[ppm_width * vi + vj] = (value > 0.0) ? log10(value) : -1.0;
    if (value > max_val) max_val = value;
    if (value < min_val && min_val > 0.0) min_val = value;
  }
  if (min_val <= 0.0) result = 1; // Matrix contains non-positive elements: log-scale corruption!
  double lmax_val = (max_val > 0.0) ? std::log10(max_val) : 0.0;
  double lmin_val = (min_val > 0.0) ? std::log10(min_val) : 0.0;

  // Write normalized pixel values.
  for (int i = 0; i < ppm_size; ++i) {
    double pix_val = (logs[i] == 1000.0) ? 0 : 255.0 * (logs[i] - lmin_val) / (lmax_val - lmin_val);
    // Avoid 0-division errors for empty matrices
    unsigned char pixel = (max_val > 0.0) ? 
      static_cast<unsigned char>(pix_val) : 0;
    ofs.put(pixel);
  }

  ofs.close();
  delete[] logs;
  return result;
}
