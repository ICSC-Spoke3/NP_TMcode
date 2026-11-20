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

/*! \file gpu_runtime.h
 *
 * \brief The few GPU runtime calls that the GPU libraries do not cover.
 *
 * MAGMA can be built on the CUDA runtime (NVIDIA GPUs) or on the HIP runtime
 * (AMD GPUs), and cuBLAS uses CUDA. The programs reach the GPU through these
 * libraries, but they also need to count the GPUs, to know their free memory
 * and to set the stack size of the GPU threads used by OpenMP offload. These
 * functions do it with the runtime of the build:
 *
 * - in MAGMA builds, the runtime MAGMA was built on (`MAGMA_HAVE_HIP` or
 *   `MAGMA_HAVE_CUDA`, defined by `magma_config.h`, which every MAGMA header
 *   includes);
 * - in cuBLAS builds without MAGMA, CUDA;
 * - otherwise, none: the functions report that no GPU runtime is available.
 *
 * No other file of the suite calls the CUDA or HIP runtime directly (the
 * cuBLAS wrapper aside).
 *
 * \ingroup support
 */

#ifndef INCLUDE_GPU_RUNTIME_H_
#define INCLUDE_GPU_RUNTIME_H_

#include <cstddef>
#include <string>

/*! \addtogroup support
 * @{
 */

/**
 * \brief Name of the GPU runtime of the build.
 *
 * \return name: `const char *` "CUDA", "HIP", or "none".
 */
const char* gpu_runtime_name();

/**
 * \brief Number of GPUs visible to the process.
 *
 * \param[out] error: `string &` Description of the failure, if the runtime cannot be used.
 * \return count: `int` The number of GPUs, or -1 if the runtime failed (or the build has none).
 */
int gpu_device_count(std::string& error);

/**
 * \brief Free memory of a GPU.
 *
 * \param device: `int` Index of the device.
 * \return bytes: `long` Free memory in bytes, or -1 if unknown.
 */
long gpu_free_memory(int device);

/**
 * \brief Set the stack size of the GPU threads of the current device.
 *
 * OpenMP offload kernels with deep call chains need a larger stack than the
 * runtime default.
 *
 * \param bytes: `size_t` Stack size of each GPU thread, in bytes.
 * \return done: `bool` Whether the runtime accepted the setting.
 */
bool gpu_set_stack_limit(size_t bytes);

//! @}
#endif
