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

/*! \file cluster.cpp
 *
 * \brief Implementation of the calculation for a cluster of spheres.
 * \ingroup prog_cluster
 */
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <hdf5.h>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#ifndef MPI_VERSION
#include <mpi.h>
#endif
#endif

#ifdef USE_NVTX
#include <nvtx3/nvToolsExt.h>
#endif

#ifdef USE_MAGMA
// In clang's OpenMP offload compilation for NVIDIA GPUs, the CUDA headers
// define __noinline__ as a macro, which breaks the libstdc++ headers (GCC 12
// and later) that use __attribute__((__noinline__)). Include those first,
// and remove the macro after the CUDA headers, unless compiling CUDA code.
#include <memory>
#include <string>
#include "magma_v2.h"
#if defined(__noinline__) && !defined(__CUDACC__)
#undef __noinline__
#endif
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
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

#ifndef INCLUDE_SPH_SUBS_H_
#include "../include/sph_subs.h"
#endif

#ifndef INCLUDE_CLU_SUBS_H_
#include "../include/clu_subs.h"
#endif

#ifndef INCLUDE_TRANSITIONMATRIX_H_
#include "../include/TransitionMatrix.h"
#endif

#ifndef INCLUDE_ALGEBRAIC_H_
#include "../include/algebraic.h"
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

#ifndef INCLUDE_RESOURCES_H_
#include "../include/resources.h"
#endif

#ifndef INCLUDE_GPU_RUNTIME_H_
#include "../include/gpu_runtime.h"
#endif

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

#ifndef INCLUDE_ITERATION_DATA_H_
#include "../include/IterationData.h"
#endif

#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

#ifdef USE_MAGMA
#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif // INCLUDE_MAGMA_CALLS
#endif // USE_MAGMA

using namespace std;

/*! \brief Main calculation loop.
 *
 *  The solution of the scattering problem for different wavelengths is an
 *  embarrasingly parallel task. This function, therefore, collects all the
 *  operations that can be independently executed by different processes,
 *  after the configuration stage and the first calculation loop have been
 *  executed.
 *
 *  \param jxi488: `int` Wavelength loop index.
 *  \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
 *  \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
 *  \param sa: `ScatteringAngles *` Pointer to a `ScatteringAngles` object.
 *  \param cid: `ClusterIterationData *` Pointer to a `ClusterIterationData` object.
 *  \param oi: `ClusterOutputInfo *` Pointer to a `ClusterOutputInfo` object.
 *  \param output_path: `const string &` Path to the output directory.
 *  \param tmatrix: `TransitionMatrix *&` Pointer that receives the T-matrix of scale
 *  JWTM, when this cycle computes it (it is written by MPI process 0 after the loop).
 *  \return jer: `int` Error code: 0 on success, non-zero if a step of the
 *  calculation (HJV, DME or the matrix inversion) failed for this wavelength.
 */
int cluster_jxi488_cycle(int jxi488, ScattererConfiguration *sconf, GeometryConfiguration *gconf, ScatteringAngles *sa, ClusterIterationData *cid, ClusterOutputInfo *oi, const string& output_path, TransitionMatrix *&tmatrix);

/*! \brief C++ implementation of CLU
 *
 *  \param config_file: `const string &` Name of the configuration file.
 *  \param data_file: `const string &` Name of the input data file.
 *  \param output_path: `const string &` Name of the directory to write the output files in.
 *  \param mpidata: `const mixMPI *` Pointer to a mixMPI data structure.
 */
void cluster(const string& config_file, const string& data_file, const string& output_path, const mixMPI *mpidata) {
  chrono::time_point<chrono::high_resolution_clock> t_start = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed;
  string message;
  string timing_name;
  FILE *timing_file;
  Logger *time_logger;
  // Checks before the calculation. MPI process 0 reads the input files and
  // writes the outputs, so it checks them; the result is shared with the
  // other processes (see agree_on_errors()), which stop with it. Other
  // processes may run on nodes that do not see the same file system.
  timing_file = NULL;
  time_logger = NULL;
  int preflight_code = 0;
  string preflight_message;
  try {
    if (mpidata->rank == 0) {
      check_input_file(config_file, "scatterer configuration");
      check_input_file(data_file, "geometry configuration");
      check_output_folder(output_path);
      timing_name = output_path + "/c_timing_mpi"+ to_string(mpidata->rank) +".log";
      timing_file = fopen(timing_name.c_str(), "w");
      if (timing_file == NULL) {
	throw OutputError("cannot create the timing log \"" + timing_name + "\"");
      }
      time_logger = new Logger(LOG_DEBG, timing_file);
    }
  } catch (const NptmError &ex) {
    preflight_code = ex.exit_code();
    preflight_message = ex.what();
  }
  agree_on_errors(mpidata, preflight_code, preflight_message);
  Logger *logger = new Logger(LOG_DEBG);
  int device_count = 0;
  // State of the GPUs of this process: a problem is not fatal here, because
  // the device policy of the calculation is known only after process 0 has
  // read the input (see check_devices()).
  int device_status = DEVICE_STATUS_OK;

#ifdef USE_MAGMA
  //===========
  // Initialise MAGMA
  //===========
  const magma_int_t d_array_max_size = 32; // TEMPORARY: can become configurable parameter
  magma_device_t *device_array = new magma_device_t[d_array_max_size];
  magma_int_t num_devices = 0;
#ifdef USE_TARGET_OFFLOAD
  // The offload kernels need a larger stack than the runtime default.
  gpu_set_stack_limit(4096);
#endif // USE TARGET_OFFLOAD
  magma_getdevices(device_array, d_array_max_size, &num_devices);
  if (num_devices < 0) num_devices = 0;
  if (num_devices > d_array_max_size) num_devices = d_array_max_size;
  device_count = (int)num_devices;
  delete[] device_array;
  message = "Proc-" + to_string(mpidata->rank) + " found " + to_string(device_count) + " GPU ";
  if (device_count > 1) message += "devices.\n";
  else message += "device.\n";
  logger->log(message, LOG_DEBG);
  logger->log("Process " + to_string(mpidata->rank) + " initializes MAGMA.\n");
  magma_int_t magma_result = magma_init();
  bool magma_ready = (magma_result == MAGMA_SUCCESS);
  if (!magma_ready) {
    logger->log(
      "Process " + to_string(mpidata->rank) + ": MAGMA initialization failed (MAGMA error code "
      + to_string(magma_result) + ").\n", LOG_DEBG
    );
    device_status = DEVICE_STATUS_INIT_FAILED;
    device_count = 0;
  } else if (device_count == 0) {
    device_status = DEVICE_STATUS_NONE_FOUND;
  } else if (!magma_device_usable(mpidata->rank % device_count)) {
    // magma_getdevices() can count GPUs that cannot be used (for instance,
    // when none is visible): check the GPU of this process.
    logger->log(
      "Process " + to_string(mpidata->rank) + ": GPU " + to_string(mpidata->rank % device_count)
      + " cannot allocate memory.\n", LOG_DEBG
    );
    device_status = DEVICE_STATUS_INIT_FAILED;
    device_count = 0;
  }
  // end MAGMA initialisation
#elif defined USE_CUBLAS
  string runtime_error;
  device_count = gpu_device_count(runtime_error);
  if (device_count < 0) {
    logger->log(
      "Process " + to_string(mpidata->rank) + ": GPU runtime initialization failed ("
      + runtime_error + ").\n", LOG_DEBG
    );
    device_status = DEVICE_STATUS_INIT_FAILED;
    device_count = 0;
  } else if (device_count == 0) {
    device_status = DEVICE_STATUS_NONE_FOUND;
  }
  logger->log(
    "Proc-" + to_string(mpidata->rank) + " found " + to_string(device_count) + " "
    + gpu_runtime_name() + " devices.\n", LOG_DEBG
  );
#endif

  // Memory and GPUs of all the processes, measured before the calculation
  // allocates its data: process 0 checks the needs of the calculation.
  ResourceSurvey survey(mpidata, device_count, device_status);

  //===========================
  // the following only happens on MPI process 0
  //===========================
  if (mpidata->rank == 0) {
    // Until the configuration is sent to the other MPI processes, they wait
    // for it: a problem found here is reported to them before stopping.
    bool setup_shared = false;
    try {
#ifdef USE_NVTX
    nvtxRangePush("Set up");
#endif
    //=======================
    // Initialise sconf from configuration file
    //=======================
    logger->log("making legacy configuration...", LOG_INFO);
    ScattererConfiguration *sconf = NULL;
    sconf = ScattererConfiguration::from_dedfb(config_file);
    sconf->write_formatted(output_path + "/c_OEDFB");
    sconf->write_binary(output_path + "/c_TEDF");
    sconf->write_binary(output_path + "/c_TEDF.hd5", "HDF5");
    // end scatterer initialisation

    //========================
    // Initialise gconf from configuration files
    //========================
    GeometryConfiguration *gconf = NULL;
    gconf = GeometryConfiguration::from_legacy(data_file);
    logger->log(" done.\n", LOG_INFO);
    //end gconf initialisation

#ifdef USE_NVTX
    nvtxRangePop();
#endif
    int s_nsph = sconf->number_of_spheres;
    int nsph = gconf->number_of_spheres;
    // Sanity check on number of sphere consistency, should always be verified
    if (s_nsph == nsph) {
      message = (gconf->offload_flag) ?
	"target offload is enabled.\n" :
	"target offload is disabled.\n";
      logger->log(message);
      char virtual_line[256];
      sprintf(virtual_line, "%.5lg m.\n", sconf->get_particle_radius(gconf));
      message = "particle radius is " + (string)virtual_line;
      logger->log(message);
      if (gconf->invert_mode == RuntimeSettings::INV_MODE_LU) {
	message = "using LU factorization for inversion.\n";
	logger->log(message);
      } else if(gconf->invert_mode == RuntimeSettings::INV_MODE_GESV) {
	message = "using GESV system solver for inversion.\n";
	logger->log(message);
      } else if(gconf->invert_mode == RuntimeSettings::INV_MODE_RBT) {
	message = "using RBT for inversion.\n";
	logger->log(message);
      }
      // else if(gconf->invert_mode == RuntimeSettings::INV_MODE_SVD) {
      // 	message = "INFO: using SVD for inversion.\n";
      // 	logger->log(message);
      // }
      // Overlapping spheres test
      double tolerance = gconf->tolerance;
      if (tolerance < 0.0) {
	tolerance = 0.1 * sconf->get_min_radius();
      }
      int *overlaps = check_overlaps(sconf, gconf, tolerance);
      if (overlaps[0] != 0) {
	message = "detected limiting separation between spheres!\n";
	logger->log(message, LOG_WARN);
	for (int oi = 0; oi < overlaps[0]; oi++) {
	  int index_0 = overlaps[1 + 2 * oi];
	  int index_1 = overlaps[2 + 2 * oi];
	  double overlap = get_overlap(sconf, gconf, index_0 - 1, index_1 - 1);
	  if (overlap > 0.0) {
	    message = "sphere " + to_string(index_0) +
	      " overlaps with sphere " + to_string(index_1) +
	      "\n";
	    logger->log(message);
	    sprintf(
	      virtual_line, "overlap is %.5lg m (tolerance is %.5lg m)\n",
	      overlap, tolerance
	    );
	    message = (string)virtual_line;
	    logger->log(message);
	    if (tolerance == 0.0 || overlap > tolerance) {
	      message = "[ERROR] sphere overlap above tolerance limit!\n";
	      logger->err(message);
	      throw(UnrecognizedConfigurationException(message));
	    } // tolerance == 0 || overlap > tolerance
	  } // ovelrap > 0
	} // oi loop
      } // overlaps[0] != 0
      delete[] overlaps;
      // Memory requirements test
      int omp_wavelength_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
      {
	int threadId = omp_get_thread_num();
	if (threadId == 0) {
	  omp_wavelength_threads = omp_get_num_threads();
	}
      }
#endif //_OPENMP
      // Order of the matrix to invert: NDIT = 2 * NSPH * LI * (LI + 2)
      const np_int matrix_order = (np_int)2 * gconf->number_of_spheres * gconf->li * (gconf->li + 2);
      const long matrix_size_bytes = (long)sizeof(dcomplex) * matrix_order * matrix_order;
      const double matrix_size_gb = matrix_size_bytes / 1024.0 / 1024.0 / 1024.0;
      // Host memory of MPI process 0, which also stores the results of all
      // the scales: the full output, one iteration data instance and one
      // inversion workspace per thread, and the one-scale output of each
      // thread other than thread 0.
      const long cid_size_bytes = ClusterIterationData::get_size(gconf, sconf);
      const long output_size_bytes = ClusterOutputInfo::get_size(sconf, gconf);
      const long scale_output_bytes = ClusterOutputInfo::get_size(sconf, gconf, 1, 1);
      // With the CPU fall-back, a thread may invert its matrix on the CPU.
      long workspace_bytes = get_inversion_workspace(matrix_order, gconf->invert_mode, gconf->refine_flag);
      if (gconf->device_policy == DEVICE_POLICY_USE_CPU_FALLBACK) {
	const long cpu_workspace = get_inversion_workspace(matrix_order, gconf->invert_mode, gconf->refine_flag, false);
	if (cpu_workspace > workspace_bytes) workspace_bytes = cpu_workspace;
      }
      const long requested_ram_bytes = output_size_bytes
	+ omp_wavelength_threads * (cid_size_bytes + workspace_bytes)
	+ (omp_wavelength_threads - 1) * scale_output_bytes;
      const double requested_ram_gb = requested_ram_bytes / 1024.0 / 1024.0 / 1024.0;
      sprintf(
	virtual_line, "iteration data require %.5lg GiB of RAM per thread.\n",
	cid_size_bytes / 1024.0 / 1024.0 / 1024.0
      );
      logger->log(virtual_line);
      sprintf(
	virtual_line, "matrix inversion requires %.5lg GiB of RAM per thread.\n",
	workspace_bytes / 1024.0 / 1024.0 / 1024.0
      );
      logger->log(virtual_line);
      sprintf(
	virtual_line, "output data require %.5lg GiB of RAM.\n",
	output_size_bytes / 1024.0 / 1024.0 / 1024.0
      );
      logger->log(virtual_line);
      sprintf(
	virtual_line, "code execution needs %.5lg GiB of RAM (%d thread(s)).\n",
	requested_ram_gb, omp_wavelength_threads
      );
      logger->log(virtual_line);
      // GPUs that are missing or failed: stop, or use the CPU (see resources.h).
      check_devices(survey, gconf->device_policy, logger);
      check_inversion_mode(gconf->invert_mode);
      // Needs of all the processes, against the memory of their nodes and
      // GPUs (see resources.h). With MAGMA or cuBLAS, each thread inverts its
      // matrix on the GPU of its process (get_device_inversion_size(): the
      // matrix and the workspaces of the inversion mode and refinement).
      // The offload pipeline builds the matrix there and inverts it in place.
      MemoryDemand demand;
      demand.full_output = output_size_bytes;
      demand.scale_output = scale_output_bytes;
      demand.per_thread = cid_size_bytes + workspace_bytes;
      demand.device_per_thread = (device_count > 0 || gconf->offload_flag) ?
	get_device_inversion_size(matrix_order, gconf->invert_mode, gconf->refine_flag) : 0;
      check_memory(survey, demand, gconf->host_ram_gb, gconf->gpu_ram_gb, gconf->memory_policy, logger);
      ScatteringAngles *p_scattering_angles = new ScatteringAngles(gconf);
      double wp = sconf->wp;
      // ClusterOutputInfo : Thread 0 of MPI process 0 allocates the memory to
      // store all the results.
      ClusterOutputInfo *p_output = new ClusterOutputInfo(sconf, gconf, mpidata);
      // Create and initialise pristine cid for MPI proc 0 and thread 0
      ClusterIterationData *cid = new ClusterIterationData(gconf, sconf, mpidata, device_count);
      const int ndi = cid->c1->nsph * cid->c1->nlim;
      np_int ndit = 2 * ndi;
      sprintf(virtual_line, "Size of matrices to invert: %ld x %ld (%.5lg GiB)\n", (long)ndit, (long)ndit, matrix_size_gb);
      message = string(virtual_line);
      logger->log(message);
      time_logger->log(message);

      str(sconf->_rcf, cid->c1);
      thdps(cid->c1->lm, cid->zpv);
      double exdc = sconf->exdc;
      double exri = sqrt(exdc);

      // Create empty virtual binary file
#ifdef USE_MAGMA
      logger->log("using MAGMA calls.\n", LOG_INFO);
#elif defined USE_CUBLAS
      logger->log("using CUBLAS calls.\n", LOG_INFO);
#elif defined USE_LAPACK
      logger->log("using LAPACK calls.\n", LOG_INFO);
#else
      logger->log("using fall-back lucin() calls.\n", LOG_INFO);
#endif
      int iavm = gconf->iavm;
      int isam = gconf->isam;
      int inpol = gconf->in_pol;
      int nxi = sconf->number_of_scales;
      int nth = p_scattering_angles->nth;
      int nths = p_scattering_angles->nths;
      int nph = p_scattering_angles->nph;
      int nphs = p_scattering_angles->nphs;

      if (sconf->idfc < 0) {
	cid->vk = cid->xip * cid->wn;
	p_output->vec_vk[0] = cid->vk;
      }

      int jxi488;
      int jer = 0;

      share_setup_status(mpidata, 0);

      setup_shared = true;

      // here go the calls that send data to be duplicated on other MPI processes from process 0 to others, using MPI broadcasts, but only if MPI is actually used
#ifdef MPI_VERSION
      if (mpidata->mpirunning) {
	gconf->mpibcast(mpidata);
	sconf->mpibcast(mpidata);	    
	cid->mpibcast(mpidata);
	p_scattering_angles->mpibcast(mpidata);
	int send_barrier = MPI_Barrier(MPI_COMM_WORLD);
      }
#endif
      // Create this variable and initialise it with a default here, so that it is defined anyway, with or without OpenMP support enabled
      int ompnumthreads = 1;
      // this is for MPI process 0 (or even if we are not using MPI at all)
      int myjxi488startoffset = 0;
      int myMPIstride = ompnumthreads;
      int myMPIblock = ompnumthreads;
      // Define here shared arrays of virtual ascii and binary files, so that thread 0 will be able to access them all later
      ClusterOutputInfo **p_outarray = NULL;
      // Errors of the threads in the loop on scales (see ScaleLoopErrors).
      ScaleLoopErrors loop_errors;
      // T-matrix of scale JWTM, if a thread of this process computes it.
      TransitionMatrix *tmatrix = NULL;
      int loop_status = 0;

#ifdef USE_NVTX
      nvtxRangePush("Parallel loop");
#endif

      //===========================================
      // open the OpenMP parallel context, so each thread can initialise its stuff
      //===========================================
#pragma omp parallel
      {
	// Create and initialise this variable here, so that if OpenMP is enabled it is local to the thread, and if OpenMP is not enabled it has a well-defiled value anyway
	int myompthread = 0;

#ifdef _OPENMP
	// If OpenMP is enabled, give actual values to myompthread and ompnumthreads, and open thread-local output files
	myompthread = omp_get_thread_num();
	if (myompthread == 0) ompnumthreads = omp_get_num_threads();
#endif

	if (myompthread == 0) {
	  // Initialise some shared variables only on thread 0
	  p_outarray = new ClusterOutputInfo*[ompnumthreads]();
	  myMPIblock = ompnumthreads;
	  myMPIstride = myMPIblock;
	}

#ifdef MPI_VERSION
	if (myompthread == 0) {
	  if (mpidata->mpirunning) {
	    // only go through this if MPI has been actually used
	    for (int rr=1; rr<mpidata->nprocs; rr++) {
	      // individually send their respective starting points to other MPI processes: they start immediately after the frequencies computed by previous processes so far
	      int remotejxi488startoffset = myMPIstride;
	      MPI_Send(&remotejxi488startoffset, 1, MPI_INT, rr, 3, MPI_COMM_WORLD);
	      int remoteMPIblock;
	      MPI_Recv(&remoteMPIblock, 1, MPI_INT, rr, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      // update myMPIstride to include the ones due to MPI process rr
	      myMPIstride += remoteMPIblock;
	    }
	    // now I know the total myMPIstride, I can send it to all processes
	    MPI_Bcast(&myMPIstride, 1, MPI_INT, 0, MPI_COMM_WORLD);
	  }
	}
#endif
	// add an omp barrier to make sure that the global variables defined by thread 0 are known to all threads below this
#pragma omp barrier

	// To test parallelism, I will now start feeding this function with "clean" copies of the parameters, so that they will not be changed by previous iterations, and each one will behave as the first one. Define all (empty) variables here, so they have the correct scope, then they get different definitions depending on thread number
	ClusterIterationData *cid_2 = NULL;
	ClusterOutputInfo *p_output_2 = NULL;
	// for threads other than the 0, create distinct copies of all relevant data, while for thread 0 just define new references / pointers to the original ones
	if (myompthread == 0) {
	  cid_2 = cid;
	  // OMP thread 0 of MPI process 0 holds the pointer to the full output structure
	  p_output_2 = p_output;
	  p_outarray[0] = p_output_2;
	} else {
	  // this is not thread 0, so do create fresh copies of all local variables
	  try {
	    cid_2 = new ClusterIterationData(*cid);
	  } catch (...) {
	    loop_errors.capture();
	  }
	}
	// make sure all threads align here: I don't want the following loop to accidentally start for thread 0, possibly modifying some variables before they are copied by all other threads
	if (myompthread==0) {
	  logger->log("Syncing OpenMP threads and starting the loop on wavelengths\n");
	  // Thread 0 of process 0 has already allocated all necessary output memory
	}
#pragma omp barrier
	// ok, now I can actually start the parallel calculations
	// we now start the loop from 1, so that the first iteration is run exactly as the others
	for (int ixi488=1; ixi488<=cid->number_of_scales; ixi488 +=myMPIstride) {
	  // the parallel loop over MPI processes covers a different set of indices for each thread
#pragma omp barrier
	  int myjxi488 = ixi488+myompthread;
	  // each thread opens new virtual files and stores their pointers in the shared array
#pragma omp barrier

	  // each MPI process handles a number of contiguous scales corresponding to its number of OMP threads at this omp level of parallelism
	  // A failing thread records its error, and all the threads stop after
	  // this block of scales (see ScaleLoopErrors).
	  if (cid_2 != NULL) {
	    try {
	      if (myjxi488 <= cid->number_of_scales) {
		if (myompthread > 0) {
		  // UPDATE: non-0 threads need to allocate memory for one scale at a time.
		  p_output_2 = new ClusterOutputInfo(sconf, gconf, mpidata, myjxi488, 1);
		  p_outarray[myompthread] = p_output_2;
		}
		int jer = cluster_jxi488_cycle(myjxi488, sconf, gconf, p_scattering_angles, cid_2, p_output_2, output_path, tmatrix);
	      } else {
		if (myompthread > 0) {
		  // If there is no input for this thread, mark to skip.
		  p_outarray[myompthread] = new ClusterOutputInfo(1);
		}
	      }
	    } catch (...) {
	      loop_errors.capture();
	    }
	  }
#pragma omp barrier
	  if (myompthread == 0) loop_status = loop_errors.agree(mpidata);
#pragma omp barrier
	  if (loop_status != 0) break;

#ifdef USE_NVTX
	  nvtxRangePush("Output concatenation");
#endif
#pragma omp barrier
	  // threads different from 0 append their virtual files to the one of thread 0, and delete them
	  if (myompthread == 0) {
	    for (int ti=1; ti<ompnumthreads; ti++) {
	      p_outarray[0]->insert(*(p_outarray[ti]));
	      delete p_outarray[ti];
	      p_outarray[ti] = NULL;
	    }
	  }
#pragma omp barrier
	  //==============================================
	  // Collect all virtual files on thread 0 of MPI process 0, and append them to disk
	  //==============================================
	  if (myompthread == 0) {
	    // thread 0 writes its virtual files, now including contributions from all threads, to disk, and deletes them
	    // p_outarray[0]->append_to_disk(output_path + "/c_OCLU");
	    // delete p_outarray[0];
	    // ******************************************************
	    // if this is the very first time, we should actually use
	    // ->write_to_disk, not ->append_to_disk
	    // ******************************************************

#ifdef MPI_VERSION
	    if (mpidata->mpirunning) {
	      // only go through this if MPI has been actually used
	      for (int rr=1; rr<mpidata->nprocs; rr++) {
		// get the data from process rr by receiving it in total memory structure
		p_outarray[0]->mpireceive(mpidata, rr);
		// get the data from process rr, creating a new virtual ascii file
		// VirtualAsciiFile *p_output = new VirtualAsciiFile(mpidata, rr);
		// append to disk and delete virtual ascii file
		// p_output->append_to_disk(output_path + "/c_OCLU");
		// delete p_output;
		
		// get the data from process rr, creating a new virtual binary file
		int test = MPI_Barrier(MPI_COMM_WORLD);
	      }
	    }
#endif
	  }
	  // end block writing to disk
#ifdef USE_NVTX
	  nvtxRangePop();
#endif
#pragma omp barrier

	} // close strided loop running on MPI processes, ixi488 loop
	// delete the shared arrays I used to make available to thread 0 the virtual files of other threads
#pragma omp barrier
	if (myompthread == 0) {
	  // After an error, the outputs of the interrupted block are left.
	  for (int ti = 1; ti < ompnumthreads; ti++) delete p_outarray[ti];
	  delete[] p_outarray;
	}
	{
	  string message = "Closing thread-local output files of thread " + to_string(myompthread) + " and syncing threads.\n";
	  logger->log(message);
	}
#ifdef USE_NVTX
	nvtxRangePop();
#endif
	delete cid_2;
      }
      delete p_scattering_angles;
      // Stop here if a thread of any process failed: no partial output files.
      loop_errors.rethrow(loop_status);
      // Only MPI process 0 writes the T-matrix (all the processes take part).
      TransitionMatrix::write_from_main_process(mpidata, tmatrix, output_path);
      delete tmatrix;
      check_output_file(output_path + "/c_OCLU.hd5", "output file");
      check_output_file(output_path + "/c_OCLU", "output file");
      p_output->write(output_path + "/c_OCLU.hd5", "HDF5");
      p_output->write(output_path + "/c_OCLU", "LEGACY");
      delete p_output;
    } // closes s_nsph == nsph check

    else { // NSPH mismatch between geometry and scatterer configurations.
      throw UnrecognizedConfigurationException(
	"Inconsistent geometry and scatterer configurations."
      );
    }
      
    delete sconf;
    delete gconf;
#ifdef USE_MAGMA
    if (magma_ready) {
      logger->log("Process " + to_string(mpidata->rank) + " finalizes MAGMA.\n");
      magma_finalize();
    }
#endif
    chrono::time_point<chrono::high_resolution_clock> t_end = chrono::high_resolution_clock::now();
    elapsed = t_end - t_start;
    string message = "Calculation lasted " + to_string(elapsed.count()) + "s.\n";
    logger->log(message);
    logger->log("Finished: output written to " + output_path + "/c_OCLU\n");
    time_logger->log(message);
    fclose(timing_file);
    delete time_logger;
    } catch (NptmError &ex) {
      if (!setup_shared) {
	share_setup_status(mpidata, ex.exit_code());
	ex.set_all_processes();
      }
      throw;
    } catch (const std::bad_alloc &ex) {
      if (!setup_shared) {
	share_setup_status(mpidata, NPTM_EXIT_MEMORY);
	MemoryError error(
	  "a memory allocation failed while preparing the calculation: reduce the model "
	  "size or the number of OpenMP threads, or use a node with more memory"
	);
	error.set_all_processes();
	throw error;
      }
      throw;
    } catch (const std::exception &ex) {
      if (!setup_shared) share_setup_status(mpidata, NPTM_EXIT_INTERNAL);
      throw;
    }
  } // end of instruction block for MPI process 0
  
    //===============================
    // instruction block for MPI processes different from 0
    //===============================
#ifdef MPI_VERSION
  else {
    // here go the code for MPI processes other than 0
    // copy gconf, sconf, cid and p_scattering_angles from MPI process 0
    receive_setup_status(mpidata);
    GeometryConfiguration *gconf = new GeometryConfiguration(mpidata);
    ScattererConfiguration *sconf = new ScattererConfiguration(mpidata);
    ClusterIterationData *cid = new ClusterIterationData(mpidata, device_count);
    ScatteringAngles *p_scattering_angles = new ScatteringAngles(mpidata);
    int recv_barrier = MPI_Barrier(MPI_COMM_WORLD);

    // Create this variable and initialise it with a default here, so that it is defined anyway, with or without OpenMP support enabled
    int ompnumthreads = 1;
    ClusterOutputInfo **p_outarray = NULL;
    // Errors of the threads in the loop on scales (see ScaleLoopErrors).
    ScaleLoopErrors loop_errors;
    // T-matrix of scale JWTM, if a thread of this process computes it.
    TransitionMatrix *tmatrix = NULL;
    int loop_status = 0;
    int myjxi488startoffset;
    int myMPIstride = ompnumthreads;
    int myMPIblock = ompnumthreads;
      
#pragma omp parallel
    {
      // Create and initialise this variable here, so that if OpenMP is enabled it is local to the thread, and if OpenMP is not enabled it has a well-defiled value anyway
      int myompthread = 0;
#ifdef _OPENMP
      // If OpenMP is enabled, give actual values to myompthread and ompnumthreads, and open thread-local output files
      myompthread = omp_get_thread_num();
      if (myompthread == 0) ompnumthreads = omp_get_num_threads();
#endif
      if (myompthread == 0) {
	// receive the start parameter from MPI process 0
	MPI_Recv(&myjxi488startoffset, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	// send my number of omp threads to process 0
	MPI_Send(&ompnumthreads, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
	// receive myMPIstride sent by MPI process 0 to all processes
	MPI_Bcast(&myMPIstride, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// allocate virtual files for each thread
	p_outarray = new ClusterOutputInfo*[ompnumthreads]();
      }
#pragma omp barrier
      // To test parallelism, I will now start feeding this function with "clean" copies of the parameters, so that they will not be changed by previous iterations, and each one will behave as the first one. Define all (empty) variables here, so they have the correct scope, then they get different definitions depending on thread number
      ClusterIterationData *cid_2 = NULL;
      ClusterOutputInfo *p_output_2 = NULL;
      // PLACEHOLDER
      // for threads other than the 0, create distinct copies of all relevant data, while for thread 0 just define new references / pointers to the original ones
      if (myompthread == 0) {
	cid_2 = cid;
      } else {
	// this is not thread 0, so do create fresh copies of all local variables
	try {
	  cid_2 = new ClusterIterationData(*cid);
	} catch (...) {
	  loop_errors.capture();
	}
      }
      // make sure all threads align here: I don't want the following loop to accidentally start for thread 0, possibly modifying some variables before they are copied by all other threads
#pragma omp barrier
      // ok, now I can actually start the parallel calculations
      // we now start the loop from 1, so that the first iteration is treated exactly on equal footing as others
      for (int ixi488=1; ixi488<=cid->number_of_scales; ixi488 +=myMPIstride) {
	// the parallel loop over MPI processes covers a different set of indices for each thread
#pragma omp barrier
	int myjxi488 = ixi488 + myjxi488startoffset + myompthread;
	// each thread opens new virtual files and stores their pointers in the shared array
#pragma omp barrier
	if (myompthread==0) logger->log("Syncing OpenMP threads and starting one iteration block on wavelengths\n");
	// ok, now I can actually start the parallel calculations
	// each MPI process handles a number of contiguous scales corresponding to its number of OMP threads at this omp level of parallelism
	// A failing thread records its error, and all the threads stop after
	// this block of scales (see ScaleLoopErrors).
	if (cid_2 != NULL) {
	  try {
	    if (myjxi488 <= cid->number_of_scales) {
	      if (myompthread > 0) {
		// UPDATE: non-0 threads need to allocate memory for one scale at a time.
		p_output_2 = new ClusterOutputInfo(sconf, gconf, mpidata, myjxi488, 1);
		p_outarray[myompthread] = p_output_2;
	      } else {
		// Thread 0 of non-zero MPI processes needs to allocate memory for the
		// output of all threads _doing something_.
		int iterstodo = cid->number_of_scales - myjxi488 + 1;
		if (iterstodo > ompnumthreads) iterstodo = ompnumthreads;
		p_output_2 = new ClusterOutputInfo(sconf, gconf, mpidata, myjxi488, iterstodo);
		p_outarray[0] = p_output_2;
	      }
	      int jer = cluster_jxi488_cycle(myjxi488, sconf, gconf, p_scattering_angles, cid_2, p_output_2, output_path, tmatrix);
	    } else {
	      p_outarray[myompthread] = new ClusterOutputInfo(1);
	    }
	  } catch (...) {
	    loop_errors.capture();
	  }
	}
#pragma omp barrier
	if (myompthread == 0) loop_status = loop_errors.agree(mpidata);
#pragma omp barrier
	if (loop_status != 0) break;
	// threads different from 0 append their virtual files to the one of thread 0, and delete them
	if (myompthread == 0) {
	  for (int ti=1; ti<ompnumthreads; ti++) {
	    p_outarray[0]->insert(*(p_outarray[ti]));
	    delete p_outarray[ti];
	    p_outarray[ti] = NULL;
	  }
	  // thread 0 sends the collected virtualfiles to thread 0 of MPI process 0, then deletes them
	  for (int rr=1; rr<mpidata->nprocs; rr++) {
	    if (rr == mpidata->rank) {
	      p_outarray[0]->mpisend(mpidata);
	      delete p_outarray[0];
	      p_outarray[0] = NULL;
	    }
	    int test = MPI_Barrier(MPI_COMM_WORLD);
	  }
	}
      } // ixi488: close strided loop running on MPI processes
      
      // Clean memory
#pragma omp barrier
      if (myompthread == 0) {
	// After an error, the outputs of the interrupted block are left.
	for (int ti = 0; ti < ompnumthreads; ti++) delete p_outarray[ti];
	delete[] p_outarray;
      }
      delete cid_2;

    } // close pragma omp parallel
    delete p_scattering_angles;
    delete sconf;
    delete gconf;
    // Stop here if a thread of any process failed.
    loop_errors.rethrow(loop_status);
    // Only MPI process 0 writes the T-matrix (all the processes take part).
    TransitionMatrix::write_from_main_process(mpidata, tmatrix, output_path);
    delete tmatrix;
#endif
#ifdef USE_MAGMA
    if (magma_ready) {
      logger->log("Process " + to_string(mpidata->rank) + " finalizes MAGMA.\n");
      magma_finalize();
    }
#endif
#ifdef MPI_VERSION
  }
#endif
  delete logger;
}

int cluster_jxi488_cycle(
  int jxi488, ScattererConfiguration *sconf, GeometryConfiguration *gconf,
  ScatteringAngles *sa, ClusterIterationData *cid, ClusterOutputInfo *output,
  const string& output_path, TransitionMatrix *&tmatrix
) {
  int nxi = sconf->number_of_scales;
  const dcomplex cc0 = 0.0 + I * 0.0;
  const double pi = acos(-1.0);
  char virtual_line[256];
  string message = "running scale iteration " + to_string(jxi488) + " of " + to_string(nxi) + ".\n";
  Logger *logger = new Logger(LOG_DEBG);
  logger->log(message);
  RuntimeSettings rs(gconf, logger);
  chrono::duration<double> elapsed;
  chrono::time_point<chrono::high_resolution_clock> interval_start, interval_end;
  int jer = 0;
  int lcalc = 0;
  int jaw = 1;
  // this is now part of cid, so don't mess with it here, just copy it by reference
  bool &is_first_scale = cid->is_first_scale;
  const int nsph = sconf->number_of_spheres;
  const np_int mxndm = gconf->mxndm;
  const int iavm = gconf->iavm;
  const int inpol = gconf->in_pol;
  const int npnt = gconf->npnt;
  const int npntts = gconf->npntts;
  const int isam = gconf->isam;
  const int jwtm = gconf->jwtm;
  int isq, ibf;
  int last_configuration;
  const int num_configs = sconf->configurations;
  const int ndirs = sa->nkks;
  int oindex = -1;
  int jindex = jxi488 - output->first_xi + 1;

#ifdef USE_NVTX
  nvtxRangePush("Prepare matrix calculation");
#endif
  double xi = sconf->get_scale(jxi488 - 1);
  double exdc = sconf->exdc;
  double exri = sqrt(exdc);
  const int idfc = (int)sconf->idfc;
  double vkarg = 0.0;
  if (idfc >= 0) {
    cid->vk = xi * cid->wn;
    vkarg = cid->vk;
    output->vec_vk[jindex - 1] = cid->vk;
    output->vec_xi[jindex - 1] = xi;
  } else {
    vkarg = xi * cid->vk;
    cid->sqsfi = 1.0 / (xi * xi);
    output->vec_vk[jindex - 1] = cid->vk;
    output->vec_xi[jindex - 1] = xi;
  }
  // Dynamic order check
  const int max_li = gconf->li;
  const int max_le = gconf->le;
  double xi_factor = (sconf->reference_variable_name.compare("XIV") == 0) ? xi : 1.0;
  const double alamb = 2.0 * pi * xi_factor / cid->vk;
  double size_par_li = 2.0 * pi * sqrt(exdc) * sconf->get_max_radius() / alamb;
  int recommended_li = 2 + (int)ceil(size_par_li + 4.05 * pow(size_par_li, 1.0 / 3.0));
  double size_par_le = 2.0 * pi * sqrt(exdc) * sconf->get_particle_radius(gconf) / alamb;
  int recommended_le = 1 + (int)ceil(size_par_le + 11.0 * pow(size_par_le, 1.0 / 3.0));
  if (recommended_li != cid->c1->li || recommended_le != cid->c1->le) {
    if (recommended_li < cid->c1->li) {
      if (gconf->dyn_order_flag > 0) {
	message = "lowering internal order from " + to_string(cid->c1->li) + " to "
	  + to_string(recommended_li) + " for scale iteration " + to_string(jxi488) + ".\n";
	logger->log(message, LOG_INFO);
      } else {
	message = "internal order " + to_string(cid->c1->li) + " too high for scale iteration "
	  + to_string(jxi488) + ".\n";
	logger->log(message, LOG_WARN);
      }
    }
    if (recommended_le < cid->c1->le) {
      if (gconf->dyn_order_flag > 0) {
	message = "lowering external order from " + to_string(cid->c1->le) + " to "
	  + to_string(recommended_le) + " for scale iteration " + to_string(jxi488) + ".\n";
	logger->log(message, LOG_INFO);
      } else {
	message = "external order " + to_string(cid->c1->le) + " too high for scale iteration "
	  + to_string(jxi488) + ".\n";
	logger->log(message, LOG_WARN);
      }
    }
    if (recommended_li < max_li || recommended_le < max_le) {
      if (gconf->dyn_order_flag > 0) {
	int new_li = (recommended_li < max_li) ? recommended_li : max_li;
	int new_le = (recommended_le < max_le) ? recommended_le : max_le;
	cid->update_orders(sconf->_rcf, new_li, new_le);
	is_first_scale = true;
	jaw = 1;
      }
    }
  }
  int li = cid->c1->li;
  int le = cid->c1->le;
  int lm = cid->c1->lm;
  np_int ndit = 2 * nsph * cid->c1->nlim;
  // End of dynamic order check
  hjv(exri, vkarg, jer, lcalc, cid->arg, cid->c1);
  if (jer != 0) {
    output->vec_ier[jindex - 1] = 1;
    output->vec_jxi[jindex - 1] = -jxi488;
    logger->log("Error in HJV for scale " + to_string(jxi488) + "!", LOG_ERRO);
    delete logger;
    return jer;
    // break; // rewrite this to go to the end of the function, to free locally allocated variables and return jer
  }
  last_configuration = 0;
  for (int i132 = 1; i132 <= nsph; i132++) {
    int iogi = cid->c1->iog[i132 - 1];
    if (iogi != i132) {
      for (int l123 = 1; l123 <= li; l123++) {
	cid->c1->rmi[l123 - 1][i132 - 1] = cid->c1->rmi[l123 - 1][iogi - 1];
	cid->c1->rei[l123 - 1][i132 - 1] = cid->c1->rei[l123 - 1][iogi - 1];
      } // l123 loop
    } else {
      last_configuration++;
      int nsh = cid->c1->nshl[last_configuration - 1];
      int ici = (nsh + 1) / 2;
      if (idfc == 0) {
	for (int ic = 0; ic < ici; ic++)
	  cid->c1->dc0[ic] = sconf->get_dielectric_constant(ic, i132 - 1, jxi488 - 1);
      } else {
	if (is_first_scale) {
	  for (int ic = 0; ic < ici; ic++)
	    cid->c1->dc0[ic] = sconf->get_dielectric_constant(ic, i132 - 1, 0);
	}
      }
      if (nsh % 2 == 0) cid->c1->dc0[ici] = exdc;
      dme(
	  cid->c1->li, i132, npnt, npntts, vkarg, exdc, exri,
	  cid->c1, jer, lcalc, cid->arg, last_configuration
	  );
      if (jer != 0) {
	output->vec_ier[jindex - 1] = 2;
	output->vec_jxi[jindex - 1] = -jxi488;
	return jer;
	//break;
      }
    }
    if (jer != 0) {
      output->vec_jxi[jindex - 1] = -jxi488;
      return jer;
      //break;
    }
  } // i132 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  if (rs.debug_am) {
    /* now, before cms, output am to p_outam0 */
    VirtualAsciiFile *outam0 = new VirtualAsciiFile();
    string outam0_name = output_path + "/c_AM0_JXI" + to_string(jxi488) + ".txt";
    sprintf(virtual_line, " AM matrix before CMS\n");
    outam0->append_line(virtual_line);
    sprintf(virtual_line, " %ld\n", (long)ndit);
    outam0->append_line(virtual_line);  
    sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
    outam0->append_line(virtual_line);
    write_dcomplex_matrix(outam0, cid->am, ndit, ndit);
    if (outam0->write_to_disk(outam0_name) != 0) {
      logger->log("could not write the debug matrix file \"" + outam0_name + "\".\n", LOG_WARN);
    }
    delete outam0;
  }
  // Without a GPU (proc_device < 0), the offload pipeline is not used.
  if (rs.use_offload && cid->proc_device >= 0) {
    // whenever rs.use_offload == true, USE_TARGET_OFFLOAD is defined, but we
    // have to check for the compilation flag in order to avoid the compiler
    // to walk the code, if the offload library functions are not available
#ifdef USE_TARGET_OFFLOAD
#ifdef USE_MAGMA
    const int ind3j_size = (cid->c1->lm + 1) * cid->c1->lm;
    const int nv3j = cid->c1->nv3j;
    const int nsphmo = nsph - 1;
    const int ncou = nsph * nsphmo;
    const int litpo = li + li + 1;
    const int litpos = litpo * litpo;
    const int lmtpo = cid->c1->lmtpo;
    const int lmtpos = cid->c1->lmtpos;
    const int rsize = li * nsph;
    const magma_int_t ndi = nsph * li * (li + 2);
    const magma_int_t nlem = cid->c1->nlem;
    const magma_int_t nlemt = nlem + nlem;
    magmaDoubleComplex* vec_am = (magmaDoubleComplex *)(cid->am[0]);
    magmaDoubleComplex *rmi = (magmaDoubleComplex *)(cid->c1->rmi[0]);
    magmaDoubleComplex *rei = (magmaDoubleComplex *)(cid->c1->rei[0]);
    double *rxx = cid->c1->rxx;
    double *ryy = cid->c1->ryy;
    double *rzz = cid->c1->rzz;
    int *ind3j = cid->c1->ind3j[0];
    double *v3j0 = cid->c1->v3j0;
    magmaDoubleComplex *vh = (magmaDoubleComplex *)(cid->c1->vh);
    magmaDoubleComplex *vyhj = (magmaDoubleComplex *)(cid->c1->vyhj);
    magmaDoubleComplex *vj0 = (magmaDoubleComplex *)(cid->c1->vj0);
    magmaDoubleComplex *vyj0 = (magmaDoubleComplex *)(cid->c1->vyj0);
    magmaDoubleComplex *gis_v = (magmaDoubleComplex *)(cid->c1->gis[0]);
    magmaDoubleComplex *gls_v = (magmaDoubleComplex *)(cid->c1->gls[0]);
    magmaDoubleComplex *sam_v = (magmaDoubleComplex *)(cid->c1->sam[0]);
    magmaDoubleComplex *vec_am0m = (magmaDoubleComplex *)(cid->c1->am0m[0]);
    magma_queue_t queue = NULL;
    magma_device_t device_id = (magma_device_t)(cid->proc_device);
    magma_queue_create(device_id, &queue);
    string gpu_failure; // description of a GPU failure, if any
#pragma omp target data map(to: rxx[0:nsph], ryy[0:nsph], rzz[0:nsph]) \
  map(to: ind3j[0:ind3j_size], v3j0[0:nv3j], vh[0:ncou*litpo]) \
  map(to: vyhj[0:ncou*litpos], vj0[0:nsph*lmtpo], vyj0[0:nsph*lmtpos]) \
  map(to: rmi[0:rsize], rei[0:rsize]) map(tofrom:vec_am[0:ndit*ndit]) \
  device(cid->proc_device)
    {
      // Initialize uninverted matrix on GPU
      magma_cms(
        vec_am, cid->c1, rxx, ryy, rzz, ind3j, v3j0, vh, vyhj, vj0, vyj0, rmi, rei, cid->proc_device
      );
      gpu_failure = magma_zinvert_resident(vec_am, (magma_int_t)ndit, jer, queue, cid->proc_device, rs);
      magma_queue_sync(queue);
      
      if (jer != 0) {
	sprintf(virtual_line, "matrix inversion returned code %d!\n", jer);
	message = virtual_line;
	logger->err(message);
      }
    } // end of target data region
    magma_queue_destroy(queue);
    if (!gpu_failure.empty()) {
      // The inversion failed before changing the matrix, which is now back
      // on the host: stop, or invert it on the CPU.
      if (rs.device_policy != DEVICE_POLICY_USE_CPU_FALLBACK) throw DeviceError(gpu_failure);
      logger->log(
	"GPU " + to_string(cid->proc_device) + " failed on scale " + to_string(jxi488) + " ("
	+ gpu_failure + "): inverting the matrix on the CPU (DEVICE_POLICY=USE_CPU_FALLBACK).\n",
	LOG_WARN
      );
      invert_matrix(cid->am, ndit, jer, output_path, jxi488, mxndm, -1, rs);
    }
    magma_ztm(
      vec_am, cid->c1, rxx, ryy, rzz, ind3j, v3j0, vh, vyhj, vj0, vyj0, sam_v,
      gis_v, gls_v, vec_am0m, cid->proc_device
    );
#else // NO_USE_MAGMA but USE_TARGET_OFFLOAD
    // TODO: implement full offload pipeline without MAGMA
    cms_flat(cid->am[0], cid->c1);
    invert_matrix(cid->am, ndit, jer, output_path, jxi488, mxndm, cid->proc_device, rs);
    ztm(cid->am, cid->c1);
#endif // USE_MAGMA
#else // NO_USE_TARGET_OFFLOAD
    throw NptmError("offload required, but not supported by this build"); // should never occur
#endif // USE_TARGET_OFFLOAD
  } else {
#ifdef USE_NVTX
    nvtxRangePush("Calculate inverted matrix");
#endif
    interval_start = chrono::high_resolution_clock::now();
    cms(cid->am[0], cid->c1);
    if (rs.debug_am) {
      VirtualAsciiFile *outam1 = new VirtualAsciiFile();
      string outam1_name = output_path + "/c_AM1_JXI" + to_string(jxi488) + ".txt";
      string outam1_ppm_name = output_path + "/c_AM1_JXI" + to_string(jxi488) + ".ppm";
      sprintf(virtual_line, " AM matrix after CMS before LUCIN\n");
      outam1->append_line(virtual_line);
      sprintf(virtual_line, " %ld\n", (long)ndit);
      outam1->append_line(virtual_line);  
      sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
      outam1->append_line(virtual_line);
      write_dcomplex_matrix(outam1, cid->am, ndit, ndit, " %5d %5d (%17.8lE,%17.8lE)\n", 1);
      if (outam1->write_to_disk(outam1_name) != 0) {
        logger->log("could not write the debug matrix file \"" + outam1_name + "\".\n", LOG_WARN);
      }
      delete outam1;
      write_matrix_as_ppm(cid->am[0], ndit, ndit, outam1_ppm_name);
    }
#ifdef USE_NVTX
    nvtxRangePop();
#endif // USE_NVTX
    interval_end = chrono::high_resolution_clock::now();
    elapsed = interval_end - interval_start;
    message = "matrix calculation for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
    logger->log(message);
    interval_start = chrono::high_resolution_clock::now();
#ifdef USE_NVTX
    nvtxRangePush("Invert the matrix");
#endif // USE_NVTX
    invert_matrix(cid->am, ndit, jer, output_path, jxi488, mxndm, cid->proc_device, rs);
    if (rs.debug_am) {
      VirtualAsciiFile *outam2 = new VirtualAsciiFile();
      string outam2_name = output_path + "/c_AM2_JXI" + to_string(jxi488) + ".txt";
      string outam2_ppm_name = output_path + "/c_AM2_JXI" + to_string(jxi488) + ".ppm";
      sprintf(virtual_line, " AM matrix after LUCIN before ZTM\n");
      outam2->append_line(virtual_line);
      sprintf(virtual_line, " %ld\n", (long)ndit);
      outam2->append_line(virtual_line);  
      sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
      outam2->append_line(virtual_line);
      write_dcomplex_matrix(outam2, cid->am, ndit, ndit);
      if (outam2->write_to_disk(outam2_name) != 0) {
        logger->log("could not write the debug matrix file \"" + outam2_name + "\".\n", LOG_WARN);
      }
      delete outam2;
      write_matrix_as_ppm(cid->am[0], ndit, ndit, outam2_ppm_name);
    }
#ifdef USE_NVTX
    nvtxRangePop();
#endif // USE_NVTX
    interval_end = chrono::high_resolution_clock::now();
    elapsed = interval_end - interval_start;
    message = "matrix inversion for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
    logger->log(message);
    if (jer != 0) {
      message = "matrix inversion ended with error code " + to_string(jer) + ".\n";
      logger->err(message);
      return jer;
      // break; // jxi488 loop: goes to memory clean
    }
    ztm(cid->am, cid->c1);
    if (rs.debug_am) {
      VirtualAsciiFile *outam3 = new VirtualAsciiFile();
      string outam3_name = output_path + "/c_AM3_JXI" + to_string(jxi488) + ".txt";
      string outam3_ppm_name = output_path + "/c_AM3_JXI" + to_string(jxi488) + ".ppm";
      sprintf(virtual_line, " AM0M matrix after ZTM\n");
      outam3->append_line(virtual_line);
      sprintf(virtual_line, " %ld\n", (long)ndit);
      outam3->append_line(virtual_line);  
      sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
      outam3->append_line(virtual_line);
      write_dcomplex_matrix(outam3, cid->c1->am0m, 2 * cid->c1->nlem, 2 * cid->c1->nlem);
      if (outam3->write_to_disk(outam3_name) != 0) {
        logger->log("could not write the debug matrix file \"" + outam3_name + "\".\n", LOG_WARN);
      }
      delete outam3;
      write_matrix_as_ppm(cid->c1->am0m[0], 2 * cid->c1->nlem, 2 * cid->c1->nlem, outam3_ppm_name);
    }
  }
  interval_start = chrono::high_resolution_clock::now();
#ifdef USE_NVTX
  nvtxRangePush("Average calculation");
#endif // USE_NVTX
  if (idfc >= 0) {
    if (jxi488 == jwtm) {
      // Keep the T-matrix: MPI process 0 writes it after the loop on scales.
      int nlemt = 2 * cid->c1->nlem;
      tmatrix = new TransitionMatrix(nlemt, lm, cid->vk, exri, cid->c1->am0m);
    }
  }
  // label 160
  double cs0 = 0.25 * cid->vk * cid->vk * cid->vk / acos(0.0);
  double csch = 0.0, qschu = 0.0, pschu = 0.0, s0mag = 0.0;
  dcomplex s0 = 0.0 + 0.0 * I;
  scr0(cid->vk, exri, cid->c1);
  double sqk = cid->vk * cid->vk * exdc;
  aps(cid->zpv, cid->c1->li, nsph, cid->c1, sqk, cid->gaps);
  rabas(inpol, cid->c1->li, nsph, cid->c1, cid->tqse, cid->tqspe, cid->tqss, cid->tqsps);
  last_configuration = 0;
  for (int i170 = 1; i170 <= nsph; i170++) {
    if (cid->c1->iog[i170 - 1] >= i170) {
      int i = i170 - 1;
      last_configuration++;
      oindex = (jindex - 1) * num_configs + last_configuration - 1;
      double albeds = cid->c1->sscs[i] / cid->c1->sexs[i];
      cid->c1->sqscs[i] *= cid->sqsfi;
      cid->c1->sqabs[i] *= cid->sqsfi;
      cid->c1->sqexs[i] *= cid->sqsfi;
      if (cid->c1->nshl[last_configuration - 1] != 1) {
	output->vec_sphere_ref_indices[oindex] = cc0;
	output->vec_sphere_sizes[oindex] = cid->c1->vsz[i];
      } else { // label 162
	output->vec_sphere_ref_indices[oindex] = cid->c1->vkt[i];
	output->vec_sphere_sizes[oindex] = cid->c1->vsz[i];
      }
      // label 164
      output->vec_sphere_scs[oindex] = cid->c1->sscs[i];
      output->vec_sphere_abs[oindex] = cid->c1->sabs[i];
      output->vec_sphere_exs[oindex] = cid->c1->sexs[i];
      output->vec_sphere_albs[oindex] = albeds;
      output->vec_sphere_sqscs[oindex] = cid->c1->sqscs[i];
      output->vec_sphere_sqabs[oindex] = cid->c1->sqabs[i];
      output->vec_sphere_sqexs[oindex] = cid->c1->sqexs[i];
      output->vec_fsas[oindex] = cid->c1->fsas[i];
      csch = 2.0 * cid->vk * cid->sqsfi / cid->c1->gcsv[i];
      s0 = cid->c1->fsas[i] * exri;
      qschu = imag(s0) * csch;
      pschu = real(s0) * csch;
      s0mag = cabs(s0) * cs0;
      output->vec_qschus[oindex] = qschu;
      output->vec_pschus[oindex] = pschu;
      output->vec_s0mags[oindex] = s0mag;
      double rapr = cid->c1->sexs[i] - cid->gaps[i];
      double cosav = cid->gaps[i] / cid->c1->sscs[i];
      output->vec_cosavs[oindex] = cosav;
      output->vec_raprs[oindex] = rapr;
      output->vec_tqek1[oindex] = cid->tqse[0][i];
      output->vec_tqsk1[oindex] = cid->tqss[0][i];
      output->vec_tqek2[oindex] = cid->tqse[1][i];
      output->vec_tqsk2[oindex] = cid->tqss[1][i];
    }
  } // i170 loop
  output->vec_fsat[jindex - 1] = cid->c1->tfsas;
  csch = 2.0 * cid->vk * cid->sqsfi / cid->c1->gcs;
  s0 = cid->c1->tfsas * exri;
  qschu = imag(s0) * csch;
  pschu = real(s0) * csch;
  s0mag = cabs(s0) * cs0;
  output->vec_qschut[jindex - 1] = qschu;
  output->vec_pschut[jindex - 1] = pschu;
  output->vec_s0magt[jindex - 1] = s0mag;
  pcrsm0(cid->vk, exri, inpol, cid->c1);
  apcra(cid->zpv, cid->c1->le, cid->c1->am0m, inpol, sqk, cid->gapm, cid->gappm);
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  interval_end = chrono::high_resolution_clock::now();
  elapsed = interval_end - interval_start;
  message = "average calculation for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
  logger->log(message);
  interval_start = chrono::high_resolution_clock::now();
#ifdef USE_NVTX
  nvtxRangePush("Angle loop");
#endif
  double th = sa->th;
  int done_dirs = 0;
  for (int jth486 = 1; jth486 <= sa->nth; jth486++) { // OpenMP portable?
    double ph = sa->ph;
    double cost = 0.0, sint = 0.0, cosp = 0.0, sinp = 0.0;
    for (int jph484 = 1; jph484 <= sa->nph; jph484++) {
      int jw = 0;
      if (sa->nk != 1 || is_first_scale) {
	upvmp(th, ph, 0, cost, sint, cosp, sinp, cid->u, cid->upmp, cid->unmp);
	if (isam >= 0) {
	  wmamp(
		0, cost, sint, cosp, sinp, inpol, cid->c1->le, 0,
		nsph, cid->argi, cid->u, cid->upmp, cid->unmp, cid->c1
		);
	  // label 182
	  apc(cid->zpv, cid->c1->le, cid->c1->am0m, cid->c1->w, sqk, cid->gap, cid->gapp);
	  raba(cid->c1->le, cid->c1->am0m, cid->c1->w, cid->tqce, cid->tqcpe, cid->tqcs, cid->tqcps);
	  jw = 1;
	}
      } else { // label 180, NK == 1 AND JXI488 > 1
	if (isam >= 0) {
	  // label 182
	  apc(cid->zpv, cid->c1->le, cid->c1->am0m, cid->c1->w, sqk, cid->gap, cid->gapp);
	  raba(cid->c1->le, cid->c1->am0m, cid->c1->w, cid->tqce, cid->tqcpe, cid->tqcs, cid->tqcps);
	  jw = 1;
	}
      }
      // label 184
      double thsl = sa->ths;
      double phsph = 0.0;
      for (int jths = 1; jths <= sa->nths; jths++) {
	double ths = thsl;
	int icspnv = 0;
	if (isam > 1) ths = th + sa->thsca;
	if (isam >= 1) {
	  phsph = 0.0;
	  if (ths < 0.0 || ths > 180.0) phsph = 180.0;
	  if (ths < 0.0) ths *= -1.0;
	  if (ths > 180.0) ths = 360.0 - ths;
	  if (phsph != 0.0) icspnv = 1;
	}
	// label 186
	double phs = sa->phs;
	for (int jphs = 1; jphs <= sa->nphs; jphs++) {
	  double costs = 0.0, sints = 0.0, cosps = 0.0, sinps = 0.0;
	  if (isam >= 1) {
	    phs = ph + phsph;
	    if (phs >= 360.0) phs -= 360.0;
	  }
	  // label 188: the scattered direction is computed again unless it
	  // cannot change, i.e. it is fixed (ISAM < 1) or there is only one
	  // incident direction (with ISAM >= 1 it follows the incident one).
	  // This differs from the original FORTRAN code and is under evaluation
	  // by its authors: see dev_notes/questions-for-authors.md, question 1.
	  bool goto190 = (sa->nks == 1) && (isam < 1 || sa->nth * sa->nph == 1)
	    && (!(is_first_scale) || jth486 > 1 || jph484 > 1);
	  if (!goto190) {
	    upvmp(ths, phs, icspnv, costs, sints, cosps, sinps, cid->us, cid->upsmp, cid->unsmp);
	    if (isam >= 0)
	      wmamp(
		    2, costs, sints, cosps, sinps, inpol, cid->c1->le,
		    0, nsph, cid->args, cid->us, cid->upsmp, cid->unsmp, cid->c1
		    );
	  }
	  // label 190
	  if (sa->nkks != 1 || is_first_scale) {
	    upvsp(
		  cid->u, cid->upmp, cid->unmp, cid->us, cid->upsmp, cid->unsmp, cid->up, cid->un, cid->ups, cid->uns,
		  cid->duk, isq, ibf, cid->scan, cid->cfmp, cid->sfmp, cid->cfsp, cid->sfsp
		  );
	    if (isam < 0) {
	      wmasp(
		    cost, sint, cosp, sinp, costs, sints, cosps, sinps,
		    cid->u, cid->up, cid->un, cid->us, cid->ups, cid->uns, isq, ibf, inpol, cid->c1->le,
		    0, nsph, cid->argi, cid->args, cid->c1
		    );
	    } else { // label 192
	      for (int i193 = 0; i193 < 3; i193++) {
		cid->up[i193] = cid->upmp[i193];
		cid->un[i193] = cid->unmp[i193];
		cid->ups[i193] = cid->upsmp[i193];
		cid->uns[i193] = cid->unsmp[i193];
	      }
	    }
	  }
	  // label 194
	  if (iavm == 1) crsm1(cid->vk, exri, cid->c1);
	  if (isam < 0) {
	    apc(cid->zpv, cid->c1->le, cid->c1->am0m, cid->c1->w, sqk, cid->gap, cid->gapp);
	    raba(cid->c1->le, cid->c1->am0m, cid->c1->w, cid->tqce, cid->tqcpe, cid->tqcs, cid->tqcps);
	    jw = 1;
	  }
	  // label 196
	  if (jaw != 0) {
	    jaw = 0;
	    mextc(cid->vk, exri, cid->c1->fsacm, cid->cextlr, cid->cext);
	    int jlr = 2;
	    for (int ilr210 = 1; ilr210 <= 2; ilr210++) {
	      int ipol = (ilr210 % 2 == 0) ? 1 : -1;
	      if (ilr210 == 2) jlr = 1;
	      double extsm = cid->c1->ecscm[ilr210 - 1];
	      double qextm = extsm * cid->sqsfi / cid->c1->gcs;
	      double extrm = extsm / cid->c1->ecs;
	      double scasm = cid->c1->scscm[ilr210 - 1];
	      double albdm = scasm / extsm;
	      double qscam = scasm * cid->sqsfi / cid->c1->gcs;
	      double scarm = scasm / cid->c1->scs;
	      double abssm = extsm - scasm;
	      double qabsm = abssm * cid->sqsfi / cid->c1->gcs;
	      double absrm = abssm / cid->c1->acs;
	      double acsecs = cid->c1->acs / cid->c1->ecs;
	      if (acsecs >= -1.0e-6 && acsecs <= 1.0e-6) absrm = 1.0;
	      dcomplex s0m = cid->c1->fsacm[ilr210 - 1][ilr210 - 1] * exri;
	      double qschum = imag(s0m) * csch;
	      double pschum = real(s0m) * csch;
	      double s0magm = cabs(s0m) * cs0;
	      double rfinrm = real(cid->c1->fsacm[ilr210 - 1][ilr210 - 1]) / real(cid->c1->tfsas);
	      double extcrm = imag(cid->c1->fsacm[ilr210 - 1][ilr210 - 1]) / imag(cid->c1->tfsas);
	      // label 208
	      if (ipol == -1) {
		output->vec_scc1[jindex - 1] = scasm;
		output->vec_abc1[jindex - 1] = abssm;
		output->vec_exc1[jindex - 1] = extsm;
		output->vec_albedc1[jindex - 1] = albdm;
		output->vec_qscamc1[jindex - 1] = qscam;
		output->vec_qabsmc1[jindex - 1] = qabsm;
		output->vec_qextmc1[jindex - 1] = qextm;
		output->vec_sccrt1[jindex - 1] = scarm;
		output->vec_abcrt1[jindex - 1] = absrm;
		output->vec_excrt1[jindex - 1] = extrm;
		// Following values are written after FSAC
		output->vec_qschuc1[jindex - 1] = qschum;
		output->vec_pschuc1[jindex - 1] = pschum;
		output->vec_s0magc1[jindex - 1] = s0magm;
	      } else if (ipol == 1) {
		output->vec_scc2[jindex - 1] = scasm;
		output->vec_abc2[jindex - 1] = abssm;
		output->vec_exc2[jindex - 1] = extsm;
		output->vec_albedc2[jindex - 1] = albdm;
		output->vec_qscamc2[jindex - 1] = qscam;
		output->vec_qabsmc2[jindex - 1] = qabsm;
		output->vec_qextmc2[jindex - 1] = qextm;
		output->vec_sccrt2[jindex - 1] = scarm;
		output->vec_abcrt2[jindex - 1] = absrm;
		output->vec_excrt2[jindex - 1] = extrm;
		// Following values are written after FSAC
		output->vec_qschuc2[jindex - 1] = qschum;
		output->vec_pschuc2[jindex - 1] = pschum;
		output->vec_s0magc2[jindex - 1] = s0magm;
	      } // ipol polarization switch
	      if (ilr210 == 1 && jlr == 2) {
		output->vec_fsac11[jindex - 1] = cid->c1->fsacm[0][0];
		output->vec_fsac21[jindex - 1] = cid->c1->fsacm[1][0];
	      } else if (ilr210 == 2 && jlr == 1) {
		output->vec_fsac22[jindex - 1] = cid->c1->fsacm[1][1];
		output->vec_fsac12[jindex - 1] = cid->c1->fsacm[0][1];
	      }
	      double rapr = cid->c1->ecscm[ilr210 - 1] - cid->gapm[2][ilr210 - 1];
	      double cosav = cid->gapm[2][ilr210 - 1] / cid->c1->scscm[ilr210 - 1];
	      double fz = rapr;
	      if (ilr210 == 1) {
		output->vec_cosavc1[jindex - 1] = cosav;
		output->vec_raprc1[jindex - 1] = rapr;
		output->vec_fkc1[jindex - 1] = fz;
	      } else if (ilr210 == 2) {
		output->vec_cosavc2[jindex - 1] = cosav;
		output->vec_raprc2[jindex - 1] = rapr;
		output->vec_fkc2[jindex - 1] = fz;
	      }
	    } // ilr210 loop
	  }
	  // label 212
	  output->vec_dir_scand[done_dirs] = cid->scan;
	  output->vec_dir_cfmp[done_dirs] = cid->cfmp;
	  output->vec_dir_sfmp[done_dirs] = cid->sfmp;
	  output->vec_dir_cfsp[done_dirs] = cid->cfsp;
	  output->vec_dir_sfsp[done_dirs] = cid->sfsp;
	  if (isam >= 0) {
	    output->vec_dir_un[3 * done_dirs] = cid->un[0];
	    output->vec_dir_un[3 * done_dirs + 1] = cid->un[1];
	    output->vec_dir_un[3 * done_dirs + 2] = cid->un[2];
	    output->vec_dir_uns[3 * done_dirs] = cid->uns[0];
	    output->vec_dir_uns[3 * done_dirs + 1] = cid->uns[1];
	    output->vec_dir_uns[3 * done_dirs + 2] = cid->uns[2];
	  } else { // label 214
	    output->vec_dir_un[3 * done_dirs] = cid->un[0];
	    output->vec_dir_un[3 * done_dirs + 1] = cid->un[1];
	    output->vec_dir_un[3 * done_dirs + 2] = cid->un[2];
	    output->vec_dir_uns[3 * done_dirs] = 0.0;
	    output->vec_dir_uns[3 * done_dirs + 1] = 0.0;
	    output->vec_dir_uns[3 * done_dirs + 2] = 0.0;
	  }
	  // label 224
	  scr2(cid->vk, vkarg, exri, cid->duk, cid->c1);
	  last_configuration = 0;
	  for (int i226 = 1; i226 <= nsph; i226++) {
	    if (cid->c1->iog[i226 - 1] >= i226) {
	      oindex = (jindex - 1) * num_configs * ndirs + num_configs * done_dirs + last_configuration;
	      last_configuration++;
	      output->vec_dir_sas11[oindex] = cid->c1->sas[i226 - 1][0][0];
	      output->vec_dir_sas21[oindex] = cid->c1->sas[i226 - 1][1][0];
	      output->vec_dir_sas12[oindex] = cid->c1->sas[i226 - 1][0][1];
	      output->vec_dir_sas22[oindex] = cid->c1->sas[i226 - 1][1][1];
	      for (int j225 = 0; j225 < 16; j225++) {
		cid->c1->vint[j225] = cid->c1->vints[i226 - 1][j225];
	      } // j225 loop
	      mmulc(cid->c1->vint, cid->cmullr, cid->cmul);
	      for (int i1 = 0; i1 < 4; i1++) {
		oindex = 16 * (jindex - 1) * num_configs * ndirs + 16 * num_configs * done_dirs + 16 * (last_configuration - 1) + 4 * i1;
		output->vec_dir_muls[oindex] = cid->cmul[i1][0];
		output->vec_dir_muls[oindex + 1] = cid->cmul[i1][1];
		output->vec_dir_muls[oindex + 2] = cid->cmul[i1][2];
		output->vec_dir_muls[oindex + 3] = cid->cmul[i1][3];
	      } // i1 loop
	      for (int i1 = 0; i1 < 4; i1++) {
		oindex = 16 * (jindex - 1) * num_configs * ndirs + 16 * num_configs * done_dirs + 16 * (last_configuration - 1) + 4 * i1;
		output->vec_dir_mulslr[oindex] = cid->cmullr[i1][0];
		output->vec_dir_mulslr[oindex + 1] = cid->cmullr[i1][1];
		output->vec_dir_mulslr[oindex + 2] = cid->cmullr[i1][2];
		output->vec_dir_mulslr[oindex + 3] = cid->cmullr[i1][3];
	      } // i1 loop
	    }
	  } // i226 loop
	  oindex = (jindex - 1) * ndirs + done_dirs;
	  output->vec_dir_sat11[oindex] = cid->c1->tsas[0][0];
	  output->vec_dir_sat21[oindex] = cid->c1->tsas[1][0];
	  output->vec_dir_sat12[oindex] = cid->c1->tsas[0][1];
	  output->vec_dir_sat22[oindex] = cid->c1->tsas[1][1];
	  pcros(cid->vk, exri, cid->c1);
	  mextc(cid->vk, exri, cid->c1->fsac, cid->cextlr, cid->cext);
	  mmulc(cid->c1->vint, cid->cmullr, cid->cmul);
	  if (jw != 0) {
	    jw = 0;
	  }
	  // label 254
	  int jlr = 2;
	  for (int ilr290 = 1; ilr290 <= 2; ilr290++) {
	    int ipol = (ilr290 % 2 == 0) ? 1 : -1;
	    if (ilr290 == 2) jlr = 1;
	    double extsec = cid->c1->ecsc[ilr290 - 1];
	    double qext = extsec * cid->sqsfi / cid->c1->gcs;
	    double extrat = extsec / cid->c1->ecs;
	    double scasec = cid->c1->scsc[ilr290 - 1];
	    double albedc = scasec / extsec;
	    double qsca = scasec * cid->sqsfi / cid->c1->gcs;
	    double scarat = scasec / cid->c1->scs;
	    double abssec = extsec - scasec;
	    double qabs = abssec * cid->sqsfi / cid->c1->gcs;
	    double absrat = 1.0;
	    double ratio = cid->c1->acs / cid->c1->ecs;
	    if (ratio < -1.0e-6 || ratio > 1.0e-6) absrat = abssec / cid->c1->acs;
	    s0 = cid->c1->fsac[ilr290 - 1][ilr290 - 1] * exri;
	    double qschu = imag(s0) * csch;
	    double pschu = real(s0) * csch;
	    s0mag = cabs(s0) * cs0;
	    double refinr = real(cid->c1->fsac[ilr290 - 1][ilr290 - 1]) / real(cid->c1->tfsas);
	    double extcor = imag(cid->c1->fsac[ilr290 - 1][ilr290 - 1]) / imag(cid->c1->tfsas);
	    // oindex remains unchanged
	    if (ipol == -1) {
	      output->vec_dir_scc1[oindex] = scasec;
	      output->vec_dir_abc1[oindex] = abssec;
	      output->vec_dir_exc1[oindex] = extsec;
	      output->vec_dir_albedc1[oindex] = albedc;
	      output->vec_dir_qscc1[oindex] = qsca;
	      output->vec_dir_qabc1[oindex] = qabs;
	      output->vec_dir_qexc1[oindex] = qext;
	      output->vec_dir_sccrt1[oindex] = scarat;
	      output->vec_dir_abcrt1[oindex] = absrat;
	      output->vec_dir_excrt1[oindex] = extrat;
	      output->vec_dir_fsac11[oindex] = cid->c1->fsac[0][0];
	      output->vec_dir_fsac21[oindex] = cid->c1->fsac[1][0];
	      output->vec_dir_sac11[oindex] = cid->c1->sac[0][0];
	      output->vec_dir_sac21[oindex] = cid->c1->sac[1][0];
	      output->vec_dir_qschuc1[oindex] = qschu;
	      output->vec_dir_pschuc1[oindex] = pschu;
	      output->vec_dir_s0magc1[oindex] = s0mag;
	    } else if (ipol == 1) {
	      output->vec_dir_scc2[oindex] = scasec;
	      output->vec_dir_abc2[oindex] = abssec;
	      output->vec_dir_exc2[oindex] = extsec;
	      output->vec_dir_albedc2[oindex] = albedc;
	      output->vec_dir_qscc2[oindex] = qsca;
	      output->vec_dir_qabc2[oindex] = qabs;
	      output->vec_dir_qexc2[oindex] = qext;
	      output->vec_dir_sccrt2[oindex] = scarat;
	      output->vec_dir_abcrt2[oindex] = absrat;
	      output->vec_dir_excrt2[oindex] = extrat;
	      output->vec_dir_fsac22[oindex] = cid->c1->fsac[1][1];
	      output->vec_dir_fsac12[oindex] = cid->c1->fsac[0][1];
	      output->vec_dir_sac22[oindex] = cid->c1->sac[1][1];
	      output->vec_dir_sac12[oindex] = cid->c1->sac[0][1];
	      output->vec_dir_qschuc2[oindex] = qschu;
	      output->vec_dir_pschuc2[oindex] = pschu;
	      output->vec_dir_s0magc2[oindex] = s0mag;
	    }
	    bool goto190 = isam >= 0 && (jths > 1 || jphs > 1);
	    if (!goto190) {
	      cid->gapv[0] = cid->gap[0][ilr290 - 1];
	      cid->gapv[1] = cid->gap[1][ilr290 - 1];
	      cid->gapv[2] = cid->gap[2][ilr290 - 1];
	      double extins = cid->c1->ecsc[ilr290 - 1];
	      double scatts = cid->c1->scsc[ilr290 - 1];
	      double rapr, cosav, fp, fn, fk, fx, fy, fz;
	      rftr(cid->u, cid->up, cid->un, cid->gapv, extins, scatts, rapr, cosav, fp, fn, fk, fx, fy, fz);
	      if (ipol == -1) {
		output->vec_dir_cosavc1[oindex] = cosav;
		output->vec_dir_raprc1[oindex] = rapr;
		output->vec_dir_flc1[oindex] = fp;
		output->vec_dir_frc1[oindex] = fn;
		output->vec_dir_fkc1[oindex] = fk;
		output->vec_dir_fxc1[oindex] = fx;
		output->vec_dir_fyc1[oindex] = fy;
		output->vec_dir_fzc1[oindex] = fz;
	      } else if (ipol == 1) {
		output->vec_dir_cosavc2[oindex] = cosav;
		output->vec_dir_raprc2[oindex] = rapr;
		output->vec_dir_flc2[oindex] = fp;
		output->vec_dir_frc2[oindex] = fn;
		output->vec_dir_fkc2[oindex] = fk;
		output->vec_dir_fxc2[oindex] = fx;
		output->vec_dir_fyc2[oindex] = fy;
		output->vec_dir_fzc2[oindex] = fz;
	      }
	      cid->tqev[0] = cid->tqce[ilr290 - 1][0];
	      cid->tqev[1] = cid->tqce[ilr290 - 1][1];
	      cid->tqev[2] = cid->tqce[ilr290 - 1][2];
	      cid->tqsv[0] = cid->tqcs[ilr290 - 1][0];
	      cid->tqsv[1] = cid->tqcs[ilr290 - 1][1];
	      cid->tqsv[2] = cid->tqcs[ilr290 - 1][2];
	      double tep, ten, tek, tsp, tsn, tsk;
	      tqr(cid->u, cid->up, cid->un, cid->tqev, cid->tqsv, tep, ten, tek, tsp, tsn, tsk);
	      if (ipol == -1) {
		output->vec_dir_tqelc1[oindex] = tep;
		output->vec_dir_tqerc1[oindex] = ten;
		output->vec_dir_tqekc1[oindex] = tek;
		output->vec_dir_tqslc1[oindex] = tsp;
		output->vec_dir_tqsrc1[oindex] = tsn;
		output->vec_dir_tqskc1[oindex] = tsk;
		output->vec_dir_tqexc1[oindex] = cid->tqce[0][0];
		output->vec_dir_tqeyc1[oindex] = cid->tqce[0][1];
		output->vec_dir_tqezc1[oindex] = cid->tqce[0][2];
		output->vec_dir_tqsxc1[oindex] = cid->tqcs[0][0];
		output->vec_dir_tqsyc1[oindex] = cid->tqcs[0][1];
		output->vec_dir_tqszc1[oindex] = cid->tqcs[0][2];
	      } else if (ipol == 1) {
		output->vec_dir_tqelc2[oindex] = tep;
		output->vec_dir_tqerc2[oindex] = ten;
		output->vec_dir_tqekc2[oindex] = tek;
		output->vec_dir_tqslc2[oindex] = tsp;
		output->vec_dir_tqsrc2[oindex] = tsn;
		output->vec_dir_tqskc2[oindex] = tsk;
		output->vec_dir_tqexc2[oindex] = cid->tqce[1][0];
		output->vec_dir_tqeyc2[oindex] = cid->tqce[1][1];
		output->vec_dir_tqezc2[oindex] = cid->tqce[1][2];
		output->vec_dir_tqsxc2[oindex] = cid->tqcs[1][0];
		output->vec_dir_tqsyc2[oindex] = cid->tqcs[1][1];
		output->vec_dir_tqszc2[oindex] = cid->tqcs[1][2];
	      }
	    }
	  } //ilr290 loop
	  for (int i = 0; i < 4; i++) {
	    oindex = 16 * (jindex - 1) * ndirs + 16 * done_dirs + 4 * i;
	    output->vec_dir_mulc[oindex] = cid->cmul[i][0];
	    output->vec_dir_mulc[oindex + 1] = cid->cmul[i][1];
	    output->vec_dir_mulc[oindex + 2] = cid->cmul[i][2];
	    output->vec_dir_mulc[oindex + 3] = cid->cmul[i][3];
	  }
	  for (int i = 0; i < 4; i++) {
	    oindex = 16 * (jindex - 1) * ndirs + 16 * done_dirs + 4 * i;
	    output->vec_dir_mulclr[oindex] = cid->cmullr[i][0];
	    output->vec_dir_mulclr[oindex + 1] = cid->cmullr[i][1];
	    output->vec_dir_mulclr[oindex + 2] = cid->cmullr[i][2];
	    output->vec_dir_mulclr[oindex + 3] = cid->cmullr[i][3];
	  }
	  if (iavm != 0) {
	    mmulc(cid->c1->vintm, cid->cmullr, cid->cmul);
	    // label 318
	    for (int i = 0; i < 4; i++) {
	      oindex = 16 * (jindex - 1) + 4 * i; // if IAVM fails, try adding directions
	      output->vec_dir_mulc[oindex] = cid->cmul[i][0];
	      output->vec_dir_mulc[oindex + 1] = cid->cmul[i][1];
	      output->vec_dir_mulc[oindex + 2] = cid->cmul[i][2];
	      output->vec_dir_mulc[oindex + 3] = cid->cmul[i][3];
	    }
	    for (int i = 0; i < 4; i++) {
	      oindex = 16 * (jindex - 1) + 4 * i; // if IAVM fails, try adding directions
	      output->vec_dir_mulclr[oindex] = cid->cmullr[i][0];
	      output->vec_dir_mulclr[oindex + 1] = cid->cmullr[i][1];
	      output->vec_dir_mulclr[oindex + 2] = cid->cmullr[i][2];
	      output->vec_dir_mulclr[oindex + 3] = cid->cmullr[i][3];
	    }
	  }
	  // label 420, continues jphs loop
	  if (isam < 1) phs += sa->phsstp;
	  done_dirs++;
	} // jphs loop, labeled 480
	if (isam <= 1) thsl += sa->thsstp;
      } // jths loop, labeled 482
      ph += sa->phstp;
    } // jph484 loop
    th += sa->thstp;
  } // jth486 loop
  // at the end of this iteration make sure to set is_first_scale to false, the next iteration will reinitialise _only_ if the order changes
  is_first_scale = 0;
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  if (jer == 0) output->vec_jxi[jindex - 1] = jxi488;
  interval_end = chrono::high_resolution_clock::now();
  elapsed = interval_end - interval_start;
  message = "angle loop for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
  logger->log(message);
  
  logger->log("finished scale iteration " + to_string(jxi488) + " of " + to_string(nxi) + ".\n");

  delete logger;

  return jer;
}

// >>> IMPLEMENTATION OF ClusterIterationData CLASS <<<
ClusterIterationData::ClusterIterationData(const GeometryConfiguration *gconf, const ScattererConfiguration *sconf, const mixMPI *mpidata, const int device_count) {
  c1 = new ParticleDescriptorCluster(gconf, sconf);
  const int ndi = c1->nsph * c1->nlim;
  const np_int ndit = 2 * ndi;
  gaps = new double[c1->nsph]();
  tqev = new double[3]();
  tqsv = new double[3]();
  tqse = new double*[2];
  tqspe = new dcomplex*[2];
  tqss = new double*[2];
  tqsps = new dcomplex*[2];
  tqce = new double*[2];
  tqcpe = new dcomplex*[2];
  tqcs = new double*[2];
  tqcps = new dcomplex*[2];
  for (int ti = 0; ti < 2; ti++) {
    tqse[ti] = new double[c1->nsph]();
    tqspe[ti] = new dcomplex[c1->nsph]();
    tqss[ti] = new double[c1->nsph]();
    tqsps[ti] = new dcomplex[c1->nsph]();
    tqce[ti] = new double[3]();
    tqcpe[ti] = new dcomplex[3]();
    tqcs[ti] = new double[3]();
    tqcps[ti] = new dcomplex[3]();
  }
  gapv = new double[3]();
  gapp = new dcomplex*[3];
  gappm = new dcomplex*[3];
  gap = new double*[3];
  gapm = new double*[3];
  for (int gi = 0; gi < 3; gi++) {
    gapp[gi] = new dcomplex[2]();
    gappm[gi] = new dcomplex[2]();
    gap[gi] = new double[2]();
    gapm[gi] = new double[2]();
  }
  u = new double[3]();
  us = new double[3]();
  un = new double[3]();
  uns = new double[3]();
  up = new double[3]();
  ups = new double[3]();
  unmp = new double[3]();
  unsmp = new double[3]();
  upmp = new double[3]();
  upsmp = new double[3]();
  argi = new double[1]();
  args = new double[1]();
  duk = new double[3]();
  cextlr = new double*[4];
  cext = new double*[4];
  cmullr = new double*[4];;
  cmul = new double*[4];
  for (int ci = 0; ci < 4; ci++) {
    cextlr[ci] = new double[4]();
    cext[ci] = new double[4]();
    cmullr[ci] = new double[4]();
    cmul[ci] = new double[4]();
  }
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[2];
      for (int zk = 0; zk < 2; zk++) {
	zpv[zi][zj][zk] = new double[2]();
      }
    }
  }
  am_vector = new dcomplex[ndit * ndit]();
  am = new dcomplex*[ndit];
  for (int ai = 0; ai < ndit; ai++) {
    am[ai] = (am_vector + ai * ndit);
  }
  
  arg = 0.0 + 0.0 * I;
  // These are suspect initializations
  scan = 0.0;
  cfmp = 0.0;
  sfmp = 0.0;
  cfsp = 0.0;
  sfsp = 0.0;
  // End of suspect initializations
  wn = sconf->wp / 3.0e8;
  xip = sconf->xip;
  sqsfi = 1.0;
  vk = 0.0;
  number_of_scales = sconf->number_of_scales;
  xiblock = (int) ceil(((double) (sconf->number_of_scales-1))/((double) mpidata->nprocs));
  lastxi = ((mpidata->rank+1) * xiblock)+1;
  firstxi = lastxi-xiblock+1;
  if (lastxi > sconf->number_of_scales) lastxi = sconf->number_of_scales;

#if defined(USE_MAGMA) || defined(USE_CUBLAS)
  // GPU of the process, or -1 if it has none (the matrices are inverted on the CPU).
  proc_device = (device_count > 0) ? mpidata->rank % device_count : -1;
#else
  proc_device = 0;
#endif

  // the first time the object is used, it will be a "first iteration", to ensure data structures are properly initialised.
  is_first_scale = 1;
}

ClusterIterationData::ClusterIterationData(const ClusterIterationData& rhs) {
  c1 = new ParticleDescriptorCluster(reinterpret_cast<ParticleDescriptorCluster &>(*(rhs.c1)));
  const int ndi = c1->nsph * c1->nlim;
  const np_int ndit = 2 * ndi;
  gaps = new double[c1->nsph]();
  for (int gi = 0; gi < c1->nsph; gi++) gaps[gi] = rhs.gaps[gi];
  tqev = new double[3]();
  tqsv = new double[3]();
  for (int ti = 0; ti < 3; ti++) {
    tqev[ti] = rhs.tqev[ti];
    tqsv[ti] = rhs.tqsv[ti];
  }
  tqse = new double*[2];
  tqspe = new dcomplex*[2];
  tqss = new double*[2];
  tqsps = new dcomplex*[2];
  tqce = new double*[2];
  tqcpe = new dcomplex*[2];
  tqcs = new double*[2];
  tqcps = new dcomplex*[2];
  for (int ti = 0; ti < 2; ti++) {
    tqse[ti] = new double[c1->nsph]();
    tqspe[ti] = new dcomplex[c1->nsph]();
    tqss[ti] = new double[c1->nsph]();
    tqsps[ti] = new dcomplex[c1->nsph]();
    for (int tj = 0; tj < c1->nsph; tj++) {
      tqse[ti][tj] = rhs.tqse[ti][tj];
      tqspe[ti][tj] = rhs.tqspe[ti][tj];
      tqss[ti][tj] = rhs.tqss[ti][tj];
      tqsps[ti][tj] = rhs.tqsps[ti][tj];
    }
    tqce[ti] = new double[3]();
    tqcpe[ti] = new dcomplex[3]();
    tqcs[ti] = new double[3]();
    tqcps[ti] = new dcomplex[3]();
    for (int tj = 0; tj < 3; tj++) {
      tqce[ti][tj] = rhs.tqce[ti][tj];
      tqcpe[ti][tj] = rhs.tqcpe[ti][tj];
      tqcs[ti][tj] = rhs.tqcs[ti][tj];
      tqcps[ti][tj] = rhs.tqcps[ti][tj];
    }
  }
  gapv = new double[3]();
  gapp = new dcomplex*[3];
  gappm = new dcomplex*[3];
  gap = new double*[3];
  gapm = new double*[3];
  for (int gi = 0; gi < 3; gi++) {
    gapv[gi] = rhs.gapv[gi];
    gapp[gi] = new dcomplex[2]();
    gappm[gi] = new dcomplex[2]();
    gap[gi] = new double[2]();
    gapm[gi] = new double[2]();
    for (int gj = 0; gj < 2; gj++) {
      gapp[gi][gj] = rhs.gapp[gi][gj];
      gappm[gi][gj] = rhs.gappm[gi][gj];
      gap[gi][gj] = rhs.gap[gi][gj];
      gapm[gi][gj] = rhs.gapm[gi][gj];
    }
  }
  u = new double[3]();
  us = new double[3]();
  un = new double[3]();
  uns = new double[3]();
  up = new double[3]();
  ups = new double[3]();
  unmp = new double[3]();
  unsmp = new double[3]();
  upmp = new double[3]();
  upsmp = new double[3]();
  duk = new double[3]();
  for (int ui = 0; ui < 3; ui++) {
    u[ui] = rhs.u[ui];
    us[ui] = rhs.us[ui];
    un[ui] = rhs.un[ui];
    uns[ui] = rhs.uns[ui];
    up[ui] = rhs.up[ui];
    ups[ui] = rhs.ups[ui];
    unmp[ui] = rhs.unmp[ui];
    unsmp[ui] = rhs.unsmp[ui];
    upmp[ui] = rhs.upmp[ui];
    upsmp[ui] = rhs.upsmp[ui];
    duk[ui] = rhs.duk[ui];
  }
  argi = new double[1]();
  args = new double[1]();
  argi[0] = rhs.argi[0];
  args[0] = rhs.args[0];
  cextlr = new double*[4];
  cext = new double*[4];
  cmullr = new double*[4];;
  cmul = new double*[4];
  for (int ci = 0; ci < 4; ci++) {
    cextlr[ci] = new double[4]();
    cext[ci] = new double[4]();
    cmullr[ci] = new double[4]();
    cmul[ci] = new double[4]();
    for (int cj = 0; cj < 4; cj++) {
      cextlr[ci][cj] = rhs.cextlr[ci][cj];
      cext[ci][cj] = rhs.cext[ci][cj];
      cmullr[ci][cj] = rhs.cmullr[ci][cj];
      cmul[ci][cj] = rhs.cmul[ci][cj];
    }
  }
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[2];
      for (int zk = 0; zk < 2; zk++) {
	zpv[zi][zj][zk] = new double[2]();
	zpv[zi][zj][zk][0] = rhs.zpv[zi][zj][zk][0];
	zpv[zi][zj][zk][1] = rhs.zpv[zi][zj][zk][1];
      }
    }
  }
  am_vector = new dcomplex[ndit * ndit]();
  // for (np_int ai = 0; ai < ndit * ndit; ai++) am_vector[ai] = rhs.am_vector[ai];
  am = new dcomplex*[ndit];
  for (np_int ai = 0; ai < ndit; ai++) {
    am[ai] = (am_vector + ai * ndit);
  }
  
  arg = rhs.arg;
  // These are suspect initializations
  scan = rhs.scan;
  cfmp = rhs.cfmp;
  sfmp = rhs.sfmp;
  cfsp = rhs.cfsp;
  sfsp = rhs.sfsp;
  // End of suspect initializations
  wn = rhs.wn;
  xip = rhs.xip;
  sqsfi = rhs.sqsfi;
  vk = rhs.vk;
  firstxi = rhs.firstxi;
  lastxi = rhs.lastxi;
  xiblock = rhs.xiblock;
  number_of_scales = rhs.number_of_scales;

  proc_device = rhs.proc_device;
  is_first_scale = rhs.is_first_scale;
}

#ifdef MPI_VERSION
ClusterIterationData::ClusterIterationData(const mixMPI *mpidata, const int device_count) {
  c1 = new ParticleDescriptorCluster(mpidata);
  const int ndi = c1->nsph * c1->nlim;
  const np_int ndit = 2 * ndi;
  gaps = new double[c1->nsph]();
  MPI_Bcast(gaps, c1->nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tqev = new double[3]();
  tqsv = new double[3]();
  MPI_Bcast(tqev, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tqsv, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tqse = new double*[2];
  tqspe = new dcomplex*[2];
  tqss = new double*[2];
  tqsps = new dcomplex*[2];
  tqce = new double*[2];
  tqcpe = new dcomplex*[2];
  tqcs = new double*[2];
  tqcps = new dcomplex*[2];
  for (int ti = 0; ti < 2; ti++) {
    tqse[ti] = new double[c1->nsph]();
    tqspe[ti] = new dcomplex[c1->nsph]();
    tqss[ti] = new double[c1->nsph]();
    tqsps[ti] = new dcomplex[c1->nsph]();
    MPI_Bcast(tqse[ti], c1->nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqspe[ti], c1->nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqss[ti], c1->nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqsps[ti], c1->nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    tqce[ti] = new double[3]();
    tqcpe[ti] = new dcomplex[3]();
    tqcs[ti] = new double[3]();
    tqcps[ti] = new dcomplex[3]();
    MPI_Bcast(tqce[ti], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqcpe[ti], 3, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqcs[ti], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqcps[ti], 3, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  }
  gapv = new double[3]();
  gapp = new dcomplex*[3];
  gappm = new dcomplex*[3];
  gap = new double*[3];
  gapm = new double*[3];
  MPI_Bcast(gapv, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (int gi = 0; gi < 3; gi++) {
    gapp[gi] = new dcomplex[2]();
    gappm[gi] = new dcomplex[2]();
    gap[gi] = new double[2]();
    gapm[gi] = new double[2]();
    MPI_Bcast(gapp[gi], 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(gappm[gi], 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(gap[gi], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(gapm[gi], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  u = new double[3]();
  us = new double[3]();
  un = new double[3]();
  uns = new double[3]();
  up = new double[3]();
  ups = new double[3]();
  unmp = new double[3]();
  unsmp = new double[3]();
  upmp = new double[3]();
  upsmp = new double[3]();
  duk = new double[3]();
  MPI_Bcast(u, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(us, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(un, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(uns, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(up, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ups, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(duk, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  argi = new double[1]();
  args = new double[1]();
  MPI_Bcast(argi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(args, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  cextlr = new double*[4];
  cext = new double*[4];
  cmullr = new double*[4];;
  cmul = new double*[4];
  for (int ci = 0; ci < 4; ci++) {
    cextlr[ci] = new double[4]();
    cext[ci] = new double[4]();
    cmullr[ci] = new double[4]();
    cmul[ci] = new double[4]();
    MPI_Bcast(cextlr[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cext[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cmullr[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cmul[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[2];
      for (int zk = 0; zk < 2; zk++) {
	zpv[zi][zj][zk] = new double[2]();
	MPI_Bcast(zpv[zi][zj][zk], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    }
  }
  am_vector = new dcomplex[ndit * ndit]();
  am = new dcomplex*[ndit];
  for (np_int ai = 0; ai < ndit; ai++) {
    am[ai] = (am_vector + ai * ndit);
    // MPI_Bcast(am[ai], ndit, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(&arg, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&scan, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&wn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xip, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sqsfi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vk, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xiblock, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&number_of_scales, 1, MPI_INT, 0, MPI_COMM_WORLD);
  lastxi = ((mpidata->rank+1) * xiblock)+1;
  firstxi = lastxi-xiblock+1;
  if (lastxi > number_of_scales) lastxi = number_of_scales;

#if defined(USE_MAGMA) || defined(USE_CUBLAS)
  // GPU of the process, or -1 if it has none (the matrices are inverted on the CPU).
  proc_device = (device_count > 0) ? mpidata->rank % device_count : -1;
#else
  proc_device = 0;
#endif
  MPI_Bcast(&is_first_scale, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
}

void ClusterIterationData::mpibcast(const mixMPI *mpidata) {
  c1->mpibcast(mpidata);
  const int ndi = c1->nsph * c1->nlim;
  const np_int ndit = 2 * ndi;
  MPI_Bcast(gaps, c1->nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tqev, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tqsv, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (int ti = 0; ti < 2; ti++) {
    MPI_Bcast(tqse[ti], c1->nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqspe[ti], c1->nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqss[ti], c1->nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqsps[ti], c1->nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqce[ti], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqcpe[ti], 3, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqcs[ti], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tqcps[ti], 3, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(gapv, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (int gi = 0; gi < 3; gi++) {
    MPI_Bcast(gapp[gi], 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(gappm[gi], 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(gap[gi], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(gapm[gi], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(u, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(us, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(un, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(uns, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(up, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ups, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(duk, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(argi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(args, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (int ci = 0; ci < 4; ci++) {
    MPI_Bcast(cextlr[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cext[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cmullr[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cmul[ci], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  for (int zi = 0; zi < c1->lm; zi++) {
    for (int zj = 0; zj < 3; zj++) {
      for (int zk = 0; zk < 2; zk++) {
	MPI_Bcast(zpv[zi][zj][zk], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    }
  }
  // since MPI expects an int argument for the number of elements to transfer in one go, transfer am one row at a time
  // for (int ai = 0; ai < ndit; ai++) {
  //   MPI_Bcast(am[ai], ndit, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  // }
  MPI_Bcast(&arg, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&scan, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&wn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xip, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sqsfi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vk, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xiblock, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&number_of_scales, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_first_scale, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
}
#endif

ClusterIterationData::~ClusterIterationData() {
  const int nsph = c1->nsph;
  delete[] am_vector;
  delete[] am;
  for (int zi = c1->lm - 1; zi > -1; zi--) {
    for (int zj = 2; zj > -1; zj--) {
      delete[] zpv[zi][zj][1];
      delete[] zpv[zi][zj][0];
      delete[] zpv[zi][zj];
    }
    delete[] zpv[zi];
  }
  delete[] zpv;
  delete c1;
  delete[] gaps;
  for (int ti = 1; ti > -1; ti--) {
    delete[] tqse[ti];
    delete[] tqss[ti];
    delete[] tqspe[ti];
    delete[] tqsps[ti];
    delete[] tqce[ti];
    delete[] tqcpe[ti];
    delete[] tqcs[ti];
    delete[] tqcps[ti];
  }
  delete[] tqse;
  delete[] tqss;
  delete[] tqspe;
  delete[] tqsps;
  delete[] tqce;
  delete[] tqcpe;
  delete[] tqcs;
  delete[] tqcps;
  delete[] tqev;
  delete[] tqsv;
  for (int gi = 2; gi > -1; gi--) {
    delete[] gapp[gi];
    delete[] gappm[gi];
    delete[] gap[gi];
    delete[] gapm[gi];
  }
  delete[] gapp;
  delete[] gappm;
  delete[] gap;
  delete[] gapm;
  delete[] gapv;
  delete[] u;
  delete[] us;
  delete[] un;
  delete[] uns;
  delete[] up;
  delete[] ups;
  delete[] unmp;
  delete[] unsmp;
  delete[] upmp;
  delete[] upsmp;
  delete[] argi;
  delete[] args;
  delete[] duk;
  for (int ci = 3; ci > -1; ci--) {
    delete[] cextlr[ci];
    delete[] cext[ci];
    delete[] cmullr[ci];
    delete[] cmul[ci];
  }
  delete[] cextlr;
  delete[] cext;
  delete[] cmullr;
  delete[] cmul;
}

long ClusterIterationData::get_size(const GeometryConfiguration *gconf, const ScattererConfiguration *sconf) {
  /*
  ParticleDescriptor *c1; double *gaps; double **tqse; dcomplex **tqspe; double **tqss; dcomplex **tqsps;
  double ****zpv; double **gapm; dcomplex **gappm; double *argi; double *args; double **gap; dcomplex **gapp;
  double **tqce; dcomplex **tqcpe; double **tqcs; dcomplex **tqcps; double *duk; double **cextlr; double **cext;
  double **cmullr; double **cmul; double *gapv; double *tqev; double *tqsv; double *u; double *us; double *un;
  double *uns; double *up; double *ups; double *unmp; double *unsmp; double *upmp; double *upsmp;
  dcomplex *am_vector; dcomplex **am;
  37 root pointers

  double scan; double cfmp; double sfmp; double cfsp; double sfsp; double sqsfi; double vk; double wn; double xip;
  9 double values
  
  dcomplex arg;
  1 dcomplex value

  int number_of_scales; int xiblock; int firstxi; int lastxi; int proc_device;
  5 int values
  
  bool is_first_scale;
  1 boolean value

  tqse = new double*[2]; tqss = new double*[2]; tqce = new double*[2]; tqcs = new double*[2];
  tqspe = new dcomplex*[2]; tqsps = new dcomplex*[2]; tqcpe = new dcomplex*[2]; tqcps = new dcomplex*[2];
  gapp = new dcomplex*[3]; gappm = new dcomplex*[3]; gap = new double*[3]; gapm = new double*[3];
  cextlr = new double*[4]; cext = new double*[4]; cmullr = new double*[4]; cmul = new double*[4];
  zpv[] = new double*[6 * LM]; am = new dcomplex*[ndit];
  (44 + 6 * LM + NDIT) long values

  tqse[] = new double[2 * c1->nsph](); tqss[] = new double[2 * c1->nsph](); tqce[] = new double[6]();
  tqcs[] = new double[6](); gaps = new double[c1->nsph](); tqev = new double[3](); tqsv = new double[3]();
  gap[] = new double[6](); gapm[] = new double[6](); gapv = new double[3](); u = new double[3]();
  us = new double[3](); un = new double[3](); uns = new double[3](); up = new double[3]();
  ups = new double[3](); unmp = new double[3](); unsmp = new double[3](); upmp = new double[3]();
  upsmp = new double[3](); argi = new double[1](); args = new double[1](); duk = new double[3]();
  cextlr[] = new double[16](); cext[] = new double[16](); cmullr[] = new double[16]();
  cmul[] = new double[16](); zpv[][] = new double[12 * lm];
  (132 + 5 * NSPH + 12 * LM) double values

  tqspe[] = new dcomplex[2 * c1->nsph](); tqsps[] = new dcomplex[2 * c1->nsph]();
  tqcpe[] = new dcomplex[6](); tqcps[] = new dcomplex[6](); gapp[] = new dcomplex[6]();
  gappm[] = new dcomplex[6](); am_vector = new dcomplex[ndit * ndit]();
  (24 + 4 * NSPH + NDIT * NDIT) dcomplex values

  */
  const int nsph = gconf->number_of_spheres;
  const int nlim = gconf->li * (gconf->li + 2);
  const int ndi = nsph * nlim;
  const np_int ndit = 2 * ndi;
  const int lm = (gconf->li > gconf->le) ? gconf->li : gconf->le;
  long result = sizeof(long) * 37;
  result += sizeof(double) * 9;
  result += sizeof(dcomplex);
  result += sizeof(int) * 5;
  result += sizeof(bool);
  result += sizeof(long) * (44 + 6 * lm + ndit);
  result += sizeof(double) * (132 + 5 * nsph + 12 * lm);
  result += sizeof(dcomplex) * (24 + 4 * nsph + ndit * ndit);
  result += ParticleDescriptorCluster::get_size(gconf, sconf);
  return result;
}

int ClusterIterationData::update_orders(double **rcf, int inner_order, int outer_order) {
  int result = 0;
  int old_lm = c1->lm;
  np_int old_ndit = 2 * c1->nsph * c1->nlim;
  ((ParticleDescriptorCluster *)c1)->update_orders(inner_order, outer_order);
  const int ndi = c1->ndi;
  const np_int ndit = (np_int)c1->ndit;
  for (int zi = 0; zi < old_lm; zi++) {
    for (int zj = 0; zj < 3; zj++) {
      for (int zk = 0; zk < 2; zk++) {
	delete[] zpv[zi][zj][zk];
      }
      delete[] zpv[zi][zj];
    }
    delete[] zpv[zi];
  }
  delete[] zpv;
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[2];
      for (int zk = 0; zk < 2; zk++) {
	zpv[zi][zj][zk] = new double[2]();
      }
    }
  }
  str(rcf, c1);
  thdps(c1->lm, zpv);
  delete[] am;
  delete[] am_vector;
  am_vector = new dcomplex[ndit * ndit]();
  am = new dcomplex*[ndit];
  for (int ai = 0; ai < ndit; ai++) {
    am[ai] = am_vector + ai * ndit;
  }
  return result;
}
// >>> END OF ClusterIterationData CLASS IMPLEMENTATION <<<
