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

/*! \file resources.h
 *
 * \brief Checks of the memory and of the GPUs available to a calculation.
 *
 * Before a calculation allocates its large data structures, the programs
 * compare the memory that it needs (estimated with the `get_size()` methods
 * of the data classes) with the memory that is available:
 *
 * - host memory, on each node: the memory that the operating system reports
 *   as available (`MemAvailable` on Linux), reduced to the limit of the
 *   control group of the process, if any (as set by batch systems such as
 *   SLURM), or the value of `HOST_RAM_GB`, if the input sets it;
 * - GPU memory, on each device used by MAGMA, cuBLAS or OpenMP offload: the
 *   free memory that the CUDA runtime reports (in MAGMA and cuBLAS builds),
 *   or the value of `GPU_RAM_GB`, if the input sets it.
 *
 * The MPI processes that run on the same node share its memory, and the ones
 * that use the same GPU share its memory, so their needs are added. When the
 * needs exceed the available memory, the memory policy of the calculation
 * decides: `STOP` (the default) stops the program with a memory error before
 * any large allocation, `FORCE_LARGE_RUN` writes a warning and goes on.
 *
 * The same survey records whether the GPUs of each process are usable: if a
 * process found no GPU, or could not initialize its GPU library, the device
 * policy decides (`check_devices()`): `STOP` (the default) stops the program
 * with a device error, `USE_CPU_FALLBACK` writes a warning and that process
 * inverts its matrices on the CPU.
 *
 * \ingroup support
 */

#ifndef INCLUDE_RESOURCES_H_
#define INCLUDE_RESOURCES_H_

#include <string>
#include <vector>

class mixMPI;
class Logger;

/*! \addtogroup support
 * @{
 */

/*! \brief What to do when a calculation needs more memory than is available.
 */
enum MemoryPolicy {
  MEMORY_POLICY_STOP = 0,           //!< Stop with a memory error (default).
  MEMORY_POLICY_FORCE_LARGE_RUN = 1 //!< Write a warning and run anyway.
};

/**
 * \brief Convert the value of a `MEMORY_POLICY=` option to a memory policy.
 *
 * \param value: `const string &` The value (`STOP` or `FORCE_LARGE_RUN`).
 * \return policy: `short` The memory policy.
 * \throws ParseError If the value is not a known policy.
 */
short parse_memory_policy(const std::string& value);

/*! \brief What to do when a GPU cannot be used.
 */
enum DevicePolicy {
  DEVICE_POLICY_STOP = 0,            //!< Stop with a device error (default).
  DEVICE_POLICY_USE_CPU_FALLBACK = 1 //!< Write a warning and use the CPU instead.
};

/**
 * \brief Convert the value of a `DEVICE_POLICY=` option to a device policy.
 *
 * \param value: `const string &` The value (`STOP` or `USE_CPU_FALLBACK`).
 * \return policy: `short` The device policy.
 * \throws ParseError If the value is not a known policy.
 */
short parse_device_policy(const std::string& value);

/*! \brief State of the GPUs of a process, at the start of a run.
 */
enum DeviceStatus {
  DEVICE_STATUS_OK = 0,          //!< The GPUs are ready, or the build does not use GPU libraries.
  DEVICE_STATUS_NONE_FOUND = 1,  //!< The build uses MAGMA or cuBLAS, but no GPU was found.
  DEVICE_STATUS_INIT_FAILED = 2  //!< The GPU library could not be initialized, or the GPU cannot be used.
};

/**
 * \brief Host memory available to this process.
 *
 * On Linux, it is the `MemAvailable` value of `/proc/meminfo`, reduced to
 * the free part of the memory limit of the control group of the process, if
 * it has one. On other systems, it is the number of free physical pages
 * reported by `sysconf()`, if available.
 *
 * \return bytes: `long` Available memory in bytes, or -1 if unknown.
 */
long get_available_host_memory();

/**
 * \brief Free memory of a GPU.
 *
 * It is known only in builds that use MAGMA or cuBLAS, from the CUDA runtime.
 *
 * \param device: `int` Index of the device.
 * \return bytes: `long` Free memory in bytes, or -1 if unknown.
 */
long get_available_device_memory(int device);

/*! \brief GPU memory set aside for each MPI process that uses a GPU, in bytes.
 *
 * Besides the data that `MemoryDemand::device_per_thread` counts, a process
 * uses GPU memory for the code loaded when a routine is first called and
 * for the work areas of the library handles (150-250 MiB measured with MAGMA
 * on an RTX 4060 Laptop, after the CUDA context).
 *
 * \note Mirrored in `src/scripts/nptm_estimates.py`, with
 * `NPTM_GPU_USABLE_FRACTION`: change them there too.
 */
#define NPTM_GPU_PROCESS_OVERHEAD (256L * 1024L * 1024L)

/*! \brief Fraction of the free GPU memory that a calculation may use.
 *
 * The rest is a margin for allocations that the estimates do not cover.
 * It does not apply to the amount given with `GPU_RAM_GB`.
 */
#define NPTM_GPU_USABLE_FRACTION 0.9

/*! \brief Memory needed by a calculation, by kind of data.
 *
 * MPI process 0 holds the results of all the scales; every process holds
 * one copy of the iteration data and of the inversion workspace per thread,
 * and the results of the scales that its threads are computing.
 */
struct MemoryDemand {
  //! \brief Bytes of the output of all the scales (MPI process 0 only).
  long full_output;
  //! \brief Bytes of the output of one scale.
  long scale_output;
  //! \brief Host bytes needed by each thread (iteration data and workspace).
  long per_thread;
  //! \brief GPU bytes needed by each thread (0 if the GPU is not used).
  long device_per_thread;

  /**
   * \brief Host memory needed by one MPI process.
   *
   * \param main_process: `bool` Whether the process is MPI process 0.
   * \param threads: `int` Number of OpenMP threads of the process.
   * \return bytes: `long` The host memory in bytes.
   *
   * \note Mirrored in `src/scripts/nptm_estimates.py`: change it there too.
   */
  long host_bytes(bool main_process, int threads) const;
};

/*! \brief Resources of the MPI processes of a run.
 *
 * All the processes build it at the same point, before the calculation
 * starts; MPI process 0 receives the data of all of them, to check the needs
 * of the calculation with `check_memory()`.
 */
class ResourceSurvey {
public:
  //! \brief Node of each process (the lowest rank of the processes on that node).
  std::vector<int> node;
  //! \brief OpenMP threads of each process.
  std::vector<int> threads;
  //! \brief GPU used by each process (-1 if none).
  std::vector<int> device;
  //! \brief Host memory available to each process, in bytes (-1 if unknown).
  std::vector<long> host_free;
  //! \brief Free memory of the GPU of each process, in bytes (-1 if unknown).
  std::vector<long> device_free;
  //! \brief State of the GPUs of each process (a `DeviceStatus` value).
  std::vector<int> device_status;

  /**
   * \brief Collect the resources of all the MPI processes.
   *
   * All the processes call it at the same point. Process `r` uses GPU
   * `r % device_count`, if `device_count` is not 0. In OpenMP offload
   * builds, a `device_count` of 0 is replaced by the number of offload
   * devices.
   *
   * \param mpidata: `const mixMPI *` MPI data of the process.
   * \param device_count: `int` Number of GPUs that the process can use.
   * \param status: `int` State of the GPUs of the process (a `DeviceStatus` value, optional, default is `DEVICE_STATUS_OK`).
   */
  ResourceSurvey(const mixMPI *mpidata, int device_count, int status = DEVICE_STATUS_OK);
};

/**
 * \brief Check that a calculation fits in the available memory.
 *
 * Called by MPI process 0, with the survey made at the start of the run and
 * the needs of the calculation. It adds the needs of the processes on each
 * node and on each GPU, and compares them with `host_ram_gb` and
 * `gpu_ram_gb`, if positive, or with the available memory otherwise. On
 * each GPU, the need includes `NPTM_GPU_PROCESS_OVERHEAD` for each process
 * that uses it, and the measured free memory is reduced to
 * `NPTM_GPU_USABLE_FRACTION`.
 *
 * \param survey: `const ResourceSurvey &` Resources of the processes.
 * \param demand: `const MemoryDemand &` Memory needed by the calculation.
 * \param host_ram_gb: `double` Host memory allowed on each node, in GiB (0 to use the available memory).
 * \param gpu_ram_gb: `double` Memory allowed on each GPU, in GiB (0 to use the free memory).
 * \param policy: `short` Memory policy (a `MemoryPolicy` value).
 * \param logger: `Logger *` Logger for the summary and the warnings.
 * \throws MemoryError If the calculation does not fit and the policy is `STOP`.
 */
void check_memory(
  const ResourceSurvey& survey, const MemoryDemand& demand, double host_ram_gb,
  double gpu_ram_gb, short policy, Logger *logger
);

/**
 * \brief Check that the GPUs of all the processes are usable.
 *
 * Called by MPI process 0, with the survey made at the start of the run.
 * If a process has no usable GPU, the policy decides: `STOP` throws a device
 * error, `USE_CPU_FALLBACK` writes a warning, and that process inverts its
 * matrices on the CPU.
 *
 * \param survey: `const ResourceSurvey &` Resources of the processes.
 * \param policy: `short` Device policy (a `DevicePolicy` value).
 * \param logger: `Logger *` Logger for the warnings.
 * \throws DeviceError If a process has no usable GPU and the policy is `STOP`.
 */
void check_devices(const ResourceSurvey& survey, short policy, Logger *logger);

/**
 * \brief Check that an allocation fits in the host memory, in a single process.
 *
 * \param bytes: `long` Memory needed, in bytes.
 * \param what: `const string &` Description of what needs the memory.
 * \param host_ram_gb: `double` Host memory allowed, in GiB (0 to use the available memory).
 * \param policy: `short` Memory policy (a `MemoryPolicy` value).
 * \param logger: `Logger *` Logger for the summary and the warnings.
 * \throws MemoryError If the allocation does not fit and the policy is `STOP`.
 */
void check_host_memory(long bytes, const std::string& what, double host_ram_gb, short policy, Logger *logger);

//! @}
#endif
