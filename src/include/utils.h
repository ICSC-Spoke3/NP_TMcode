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
 * \ingroup support
 */

#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <exception>
#include <string>

#ifndef INCLUDE_TYPES_H_
#include "types.h"
#endif

class mixMPI;
class VirtualAsciiFile;

/*! \addtogroup support
 * @{
 */

/**
 * \brief Estimate the host memory used by one matrix inversion.
 *
 * Returns the bytes that one call of `invert_matrix()` allocates on the host
 * for an n x n matrix, besides the matrix itself, with the inversion back end
 * selected at build time (MAGMA, cuBLAS, LAPACK or the built-in one), the
 * inversion mode and the refinement setting. Each thread that inverts a
 * matrix needs this workspace.
 *
 * In builds with MAGMA or cuBLAS, `on_gpu` selects between the GPU back end
 * and the CPU one that inverts the matrices of a process without GPU, or the
 * ones that a GPU failed to invert (`DEVICE_POLICY=USE_CPU_FALLBACK`).
 *
 * \param n: `np_int` Number of rows (and columns) of the matrix.
 * \param invert_mode: `short` Inversion mode (one of the `RuntimeSettings::INV_MODE_*` values).
 * \param refine: `bool` Whether iterative refinement is enabled.
 * \param on_gpu: `bool` Whether the GPU back end inverts the matrix (optional, default is true).
 * \return size: `long` Workspace size in bytes.
 *
 * \note Mirrored in `src/scripts/nptm_estimates.py`: change it there too.
 */
long get_inversion_workspace(np_int n, short invert_mode, bool refine, bool on_gpu = true);

/**
 * \brief Device memory needed to invert a matrix on the GPU.
 *
 * Peak GPU memory used by one call of `invert_matrix()` on an n x n matrix,
 * including the matrix itself, with the GPU back end selected at build time,
 * the inversion mode and the refinement setting:
 *
 * - MAGMA, LU: the matrix and the LU workspaces (panel and `zgetri` work,
 *   counted with the largest MAGMA block size, 512); with refinement, the
 *   matrix and the three matrices of `magma_newton()`;
 * - MAGMA, GESV: the matrix and the identity (plus the panel workspace);
 *   with refinement, three more matrices;
 * - MAGMA, RBT: about two matrices allocated by `magma_zgesv_rbt()`, five
 *   with its iterative refinement;
 * - MAGMA, SVD: the U and Vt factors and the pseudo-inverse; with
 *   refinement, three more matrices;
 * - cuBLAS: the matrix, the inverse and the cuSOLVER workspace (counted as
 *   512 columns); with refinement, two more matrices. These buffers are kept
 *   by each thread between inversions.
 *
 * In builds without MAGMA or cuBLAS, the matrices are inverted on the CPU and
 * the result is 0. The offload pipeline of `np_cluster` builds the matrix on
 * the GPU and inverts it there with the same MAGMA routines; its other
 * arrays are small and not counted.
 *
 * \param n: `np_int` Number of rows (and columns) of the matrix.
 * \param invert_mode: `short` Inversion mode (one of the `RuntimeSettings::INV_MODE_*` values).
 * \param refine: `bool` Whether iterative refinement is enabled.
 * \return size: `long` Device memory in bytes, for each thread that inverts a matrix.
 *
 * \note Mirrored in `src/scripts/nptm_estimates.py`: change it there too.
 */
long get_device_inversion_size(np_int n, short invert_mode, bool refine);

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

//! @}
/**
 * \brief Check that an input file exists and can be read.
 *
 * \param file_name: `const string &` Path of the file.
 * \param description: `const string &` What the file contains, for the message.
 * \throws InputError If the file does not exist or cannot be read.
 */
void check_input_file(const std::string& file_name, const std::string& description);

/**
 * \brief Check that an output folder exists and can be written.
 *
 * \param path: `const string &` Path of the folder.
 * \throws OutputError If the folder does not exist, is not a folder or is
 * not writable.
 */
void check_output_folder(const std::string& path);

/**
 * \brief Check that an output file can be created.
 *
 * The file is created (or emptied), as the writer that follows would do.
 *
 * \param file_name: `const string &` Path of the file.
 * \param description: `const string &` What the file contains, for the message.
 * \throws OutputError If the file cannot be created (for instance, if a
 * folder with its name exists).
 */
void check_output_file(const std::string& file_name, const std::string& description);

/**
 * \brief Make all the MPI processes stop if any of them found a problem.
 *
 * Every process calls this function at the same point, with the exit status
 * of the problem it found (0 if none). If any status is not 0, every process
 * throws: the ones that found a problem throw an exception of its category
 * with their message, the others a `PeerProcessError`. All the exceptions
 * are marked as known to all the processes, so that MPI can be finalized
 * normally. Without MPI, the function throws only if the status is not 0.
 *
 * \param mpidata: `const mixMPI *` MPI data of the process.
 * \param local_code: `int` Exit status of the problem found by this process (0 if none).
 * \param local_message: `const string &` Description of the problem.
 */
void agree_on_errors(const mixMPI *mpidata, int local_code, const std::string& local_message);

/**
 * \brief Send the result of the set-up of MPI process 0 to the other processes.
 *
 * MPI process 0 reads and checks the configuration, then sends it to the
 * other processes. It calls this function just before, with 0 if the set-up
 * succeeded or the exit status of the problem otherwise, and the other
 * processes call `receive_setup_status()` before they receive the
 * configuration. Without MPI, the function does nothing.
 *
 * \param mpidata: `const mixMPI *` MPI data of the process.
 * \param status: `int` 0 or the exit status of the problem.
 */
void share_setup_status(const mixMPI *mpidata, int status);

/**
 * \brief Receive the result of the set-up of MPI process 0.
 *
 * \param mpidata: `const mixMPI *` MPI data of the process.
 * \throws PeerProcessError If the set-up of process 0 failed.
 */
void receive_setup_status(const mixMPI *mpidata);

/*! \brief Errors of the threads and processes of a parallel loop on scales.
 *
 * Exceptions cannot leave an OpenMP parallel region, and an MPI process that
 * stops while the others wait for its results leaves them waiting forever.
 * In the loops on scales, each thread catches its exceptions and records
 * them with `capture()`. After each block of scales, thread 0 of every
 * process calls `agree()`, at the same point of all the processes, and all
 * the threads leave the loop if any process failed. After the parallel
 * region, `rethrow()` throws the error again in the process that found it,
 * and a `PeerProcessError` in the others, all marked as known to all the
 * processes, so that MPI is finalized normally.
 */
class ScaleLoopErrors {
protected:
  //! \brief The first exception captured by a thread of this process.
  std::exception_ptr _error;
  //! \brief Exit status of the captured exception (0 if none).
  int _status;

public:
  /**
   * \brief Create a record without errors.
   */
  ScaleLoopErrors();

  /**
   * \brief Record the exception that is being handled.
   *
   * Call it inside a `catch` block, from any thread: only the first
   * exception of the process is kept.
   */
  void capture() noexcept;

  /**
   * \brief Exit status of the error of this process.
   *
   * \return status: `int` The exit status, or 0 if no thread failed.
   */
  int status() const { return _status; }

  /**
   * \brief Combine the errors of all the processes.
   *
   * One thread of every MPI process calls it, at the same point.
   *
   * \param mpidata: `const mixMPI *` MPI data of the process.
   * \return status: `int` The highest exit status of the processes, or 0 if none failed.
   */
  int agree(const mixMPI *mpidata);

  /**
   * \brief Throw the error of the loop, if any.
   *
   * \param global_status: `int` The status returned by `agree()`.
   * \throws NptmError The error of this process, or a `PeerProcessError` if another process failed.
   */
  void rethrow(int global_status);
};

/**
 * \brief Report an exception that stops a program and return its exit status.
 *
 * Programs call this function from the handlers of `main()`. It logs the
 * description of the problem (not for a `PeerProcessError`, which the
 * process that found the problem reports) with the rank of the process when
 * several MPI processes run.
 *
 * \param ex: `const std::exception &` The exception.
 * \param mpidata: `const mixMPI *` MPI data of the process.
 * \return exit_code: `int` The exit status of the program.
 */
int report_stopping_error(const std::exception& ex, const mixMPI *mpidata);

#endif
