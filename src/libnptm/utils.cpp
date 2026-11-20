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
 * \ingroup support
 */

#include <exception>
#include <fstream>
#include <new>
#include <cstring>
#include <unistd.h>
#include <sys/stat.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

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

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif

using namespace std;

//! \cond
// Host workspace of a CPU inversion: LAPACK (LU or GESV, the other modes
// fall back to LU), or lucin().
static long cpu_inversion_workspace(np_int n, short invert_mode, bool refine) {
  const long matrix = (long)sizeof(dcomplex) * n * n;
  long result = 0;
#ifdef USE_LAPACK
  if (invert_mode != RuntimeSettings::INV_MODE_GESV) invert_mode = RuntimeSettings::INV_MODE_LU;
  result += sizeof(np_int) * n; // pivots
  if (invert_mode == RuntimeSettings::INV_MODE_GESV) result += matrix; // identity matrix
  if (refine) {
    // copy of the original matrix, and the product, residual and unrefined
    // inverse used by the Newton refinement
    result += 4 * matrix + (long)sizeof(dcomplex) * n;
  }
#else
  result += 2 * (long)sizeof(double) * n; // work vectors of lucin()
#endif
  return result;
}
//! \endcond

long get_inversion_workspace(np_int n, short invert_mode, bool refine, bool on_gpu) {
  // Host memory allocated by one call of invert_matrix() on an n x n matrix,
  // besides the matrix itself, for the back end selected at build time.
  const long matrix = (long)sizeof(dcomplex) * n * n;
  long result = 0;
#ifdef USE_MAGMA
  if (!on_gpu) return cpu_inversion_workspace(n, invert_mode, refine);
  result += sizeof(magma_int_t) * n; // pivots
  if (invert_mode == RuntimeSettings::INV_MODE_GESV) {
    result += matrix; // identity matrix
  } else if (invert_mode == RuntimeSettings::INV_MODE_RBT) {
    result += 2 * matrix; // identity matrix and copy of the matrix
  } else if (invert_mode == RuntimeSettings::INV_MODE_SVD) {
    // copies of the matrix and of the U and Vt factors, plus the real work
    // array of 5 n^2 + 7 n doubles
    result += 3 * matrix + (long)sizeof(double) * (5 * n * n + 7 * n);
  }
  if (refine) result += (long)sizeof(dcomplex) * n; // refinement runs on the device
#elif defined USE_CUBLAS
  if (!on_gpu) return cpu_inversion_workspace(n, invert_mode, refine);
  result += sizeof(dcomplex); // the inversion runs on the device
#else
  result = cpu_inversion_workspace(n, invert_mode, refine);
#endif
  return result;
}

long get_device_inversion_size(np_int n, short invert_mode, bool refine) {
  const long matrix = (long)sizeof(dcomplex) * n * n;
  const long vector = (long)sizeof(dcomplex) * n;
  // Panel and work arrays of the LU routines: n x nb, with the largest
  // block size that MAGMA chooses.
  const long block_columns = 512;
  long result = 0;
#ifdef USE_MAGMA
  if (invert_mode == RuntimeSettings::INV_MODE_GESV) {
    result = 2 * matrix + block_columns * vector; // matrix, identity (solution), panel
    if (refine) result = 5 * matrix + vector; // magma_newton() runs while both are allocated
  } else if (invert_mode == RuntimeSettings::INV_MODE_RBT) {
    // Matrix and right-hand sides inside magma_zgesv_rbt(), plus copies and
    // the solution when it refines.
    result = (refine ? 5 : 2) * matrix + 4 * vector;
  } else if (invert_mode == RuntimeSettings::INV_MODE_SVD) {
    result = 3 * matrix; // U, Vt and the pseudo-inverse
    if (refine) result = 6 * matrix + vector;
  } else { // LU
    result = matrix + 2 * block_columns * vector; // matrix, zgetrf panel, zgetri work
    // The work array is released before magma_newton() allocates its three matrices.
    if (refine && 4 * matrix + vector > result) result = 4 * matrix + vector;
  }
#elif defined USE_CUBLAS
  (void)invert_mode; // the cuBLAS build inverts by LU only
  result = 2 * matrix + block_columns * vector + (long)sizeof(int64_t) * n; // matrix, inverse, workspace, pivots
  if (refine) result += 2 * matrix; // residual and correction
#else
  (void)n;
  (void)invert_mode;
  (void)refine;
  (void)matrix;
  (void)vector;
  (void)block_columns;
#endif
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

void check_input_file(const std::string& file_name, const std::string& description) {
  struct stat info;
  if (stat(file_name.c_str(), &info) != 0) {
    throw InputError(description + " file \"" + file_name + "\" does not exist");
  }
  if (S_ISDIR(info.st_mode)) {
    throw InputError(description + " file \"" + file_name + "\" is a folder");
  }
  if (access(file_name.c_str(), R_OK) != 0) {
    throw InputError(description + " file \"" + file_name + "\" cannot be read (check its permissions)");
  }
}

void check_output_file(const std::string& file_name, const std::string& description) {
  std::ofstream probe(file_name, std::ios::binary | std::ios::out | std::ios::trunc);
  if (!probe.is_open()) {
    throw OutputError("cannot write the " + description + " \"" + file_name + "\"");
  }
}

void check_output_folder(const std::string& path) {
  struct stat info;
  if (stat(path.c_str(), &info) != 0) {
    throw OutputError("output folder \"" + path + "\" does not exist (create it before running)");
  }
  if (!S_ISDIR(info.st_mode)) {
    throw OutputError("output path \"" + path + "\" is not a folder");
  }
  if (access(path.c_str(), W_OK | X_OK) != 0) {
    throw OutputError("output folder \"" + path + "\" is not writable (check its permissions)");
  }
}

void agree_on_errors(const mixMPI *mpidata, int local_code, const std::string& local_message) {
  int global_code = local_code;
#ifdef MPI_VERSION
  if (mpidata->mpirunning) {
    MPI_Allreduce(&local_code, &global_code, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  if (global_code == 0) return;
  if (local_code != 0) {
    // Rethrow as the category of the problem, marked as known to all processes.
    switch (local_code) {
    case NPTM_EXIT_USAGE: { UsageError ex(local_message); ex.set_all_processes(); throw ex; }
    case NPTM_EXIT_INPUT: { InputError ex(local_message); ex.set_all_processes(); throw ex; }
    case NPTM_EXIT_PARSE: { ParseError ex(local_message); ex.set_all_processes(); throw ex; }
    case NPTM_EXIT_OUTPUT: { OutputError ex(local_message); ex.set_all_processes(); throw ex; }
    case NPTM_EXIT_MEMORY: { MemoryError ex(local_message); ex.set_all_processes(); throw ex; }
    case NPTM_EXIT_DEVICE: { DeviceError ex(local_message); ex.set_all_processes(); throw ex; }
    default: { NptmError ex(local_message, local_code); ex.set_all_processes(); throw ex; }
    }
  }
  throw PeerProcessError(global_code);
}

void share_setup_status(const mixMPI *mpidata, int status) {
#ifdef MPI_VERSION
  if (mpidata->mpirunning) {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif
}

void receive_setup_status(const mixMPI *mpidata) {
  int status = 0;
#ifdef MPI_VERSION
  if (mpidata->mpirunning) {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif
  if (status != 0) throw PeerProcessError(status);
}

int report_stopping_error(const std::exception& ex, const mixMPI *mpidata) {
  Logger logger(LOG_DEBG);
  std::string where = (mpidata != NULL && mpidata->nprocs > 1) ?
    "process " + std::to_string(mpidata->rank) + ": " : "";
  const PeerProcessError *peer = dynamic_cast<const PeerProcessError *>(&ex);
  if (peer != NULL) return peer->exit_code();
  const NptmError *nptm = dynamic_cast<const NptmError *>(&ex);
  if (nptm != NULL) {
    logger.err(where + nptm->category() + ": " + nptm->what() + "\n");
    return nptm->exit_code();
  }
  if (dynamic_cast<const std::bad_alloc *>(&ex) != NULL) {
    logger.err(
      where + "memory error: a memory allocation failed. The calculation does not fit in "
      "the available memory: reduce the model size or the number of OpenMP threads, or "
      "use a node with more memory.\n"
    );
    return NPTM_EXIT_MEMORY;
  }
  logger.err(where + "internal error: " + ex.what() + "\n");
  return NPTM_EXIT_INTERNAL;
}

ScaleLoopErrors::ScaleLoopErrors() : _error(nullptr), _status(0) { }

void ScaleLoopErrors::capture() noexcept {
  std::exception_ptr current = std::current_exception();
  int code = NPTM_EXIT_INTERNAL;
  try {
    std::rethrow_exception(current);
  } catch (const NptmError &ex) {
    code = ex.exit_code();
  } catch (const std::bad_alloc &) {
    code = NPTM_EXIT_MEMORY;
  } catch (...) {
    code = NPTM_EXIT_INTERNAL;
  }
#pragma omp critical(nptm_scale_loop_errors)
  {
    if (_status == 0) {
      _error = current;
      _status = code;
    }
  }
}

int ScaleLoopErrors::agree(const mixMPI *mpidata) {
  int global_status = _status;
#ifdef MPI_VERSION
  if (mpidata != NULL && mpidata->mpirunning) {
    MPI_Allreduce(&_status, &global_status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  return global_status;
}

void ScaleLoopErrors::rethrow(int global_status) {
  if (_error) {
    try {
      std::rethrow_exception(_error);
    } catch (NptmError &ex) {
      ex.set_all_processes();
      throw;
    } catch (const std::bad_alloc &) {
      MemoryError error(
        "a memory allocation failed during the calculation: reduce the model size or "
        "the number of OpenMP threads, or use a node with more memory"
      );
      error.set_all_processes();
      throw error;
    } catch (const std::exception &ex) {
      NptmError error(ex.what());
      error.set_all_processes();
      throw error;
    } catch (...) {
      NptmError error("unknown error during the calculation");
      error.set_all_processes();
      throw error;
    }
  }
  if (global_status != 0) throw PeerProcessError(global_status);
}
