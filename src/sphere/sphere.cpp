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

/*! \file sphere.cpp
 *
 * \brief Implementation of the single sphere calculation.
 * \ingroup prog_sphere
 */
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

#ifndef INCLUDE_TRANSITIONMATRIX_H_
#include "../include/TransitionMatrix.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

#ifndef INCLUDE_ITERATION_DATA_H_
#include "../include/IterationData.h"
#endif

#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif

#ifndef INCLUDE_RESOURCES_H_
#include "../include/resources.h"
#endif

using namespace std;

/*! \brief Main calculation loop.
 *
 *  \param jxi488: `int` Wavelength loop index.
 *  \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
 *  \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
 *  \param sa: `ScatteringAngles *` Pointer to a `ScatteringAngles` object.
 *  \param sid: `SphereIterationData *` Pointer to a `SphereIterationData` object.
 *  \param oi: `SphereOutputInfo *` Pointer to a `SphereOutputInfo` object.
 *  \param output_path: `const string &` Path to the output directory.
 *  \param tmatrix: `TransitionMatrix *&` Pointer that receives the T-matrix of scale
 *  JWTM, when this cycle computes it (it is written by MPI process 0 after the loop).
 *  \return jer: `int` Error code: 0 on success, non-zero if a step of the
 *  calculation (DME) failed for this wavelength.
 */
int sphere_jxi488_cycle(
  int jxi488, ScattererConfiguration *sconf, GeometryConfiguration *gconf,
  ScatteringAngles *sa, SphereIterationData *sid, SphereOutputInfo *oi,
  const string& output_path, TransitionMatrix *&tmatrix
);

/*! \brief C++ implementation of SPH
 *
 *  \param config_file: `string` Name of the configuration file.
 *  \param data_file: `string` Name of the input data file.
 *  \param output_path: `string` Directory to write the output files in.
 *  \param mpidata: `const mixMPI *` Pointer to a mixMPI data structure.
 */
void sphere(const string& config_file, const string& data_file, const string& output_path, const mixMPI *mpidata) {
  Logger *logger = new Logger(LOG_INFO);
  int device_count = 0;

  // Checks before the calculation. MPI process 0 reads the input files and
  // writes the outputs, so it checks them; the result is shared with the
  // other processes (see agree_on_errors()), which stop with it. Other
  // processes may run on nodes that do not see the same file system.
  int preflight_code = 0;
  string preflight_message;
  try {
    if (mpidata->rank == 0) {
      check_input_file(config_file, "scatterer configuration");
      check_input_file(data_file, "geometry configuration");
      check_output_folder(output_path);
    }
  } catch (const NptmError &ex) {
    preflight_code = ex.exit_code();
    preflight_message = ex.what();
  }
  agree_on_errors(mpidata, preflight_code, preflight_message);

  // Memory of all the processes, measured before the calculation allocates
  // its data: process 0 checks the needs of the calculation.
  ResourceSurvey survey(mpidata, device_count);

  //===========================
  // the following only happens on MPI process 0
  //===========================
  if (mpidata->rank == 0) {
    // Until the configuration is sent to the other MPI processes, they wait
    // for it: a problem found here is reported to them before stopping.
    bool setup_shared = false;
    try {
    logger->log("making legacy configuration...");
    ScattererConfiguration *sconf = NULL;
    sconf = ScattererConfiguration::from_dedfb(config_file);
    sconf->write_formatted(output_path + "/c_OEDFB");
    sconf->write_binary(output_path + "/c_TEDF");
    sconf->write_binary(output_path + "/c_TEDF.hd5", "HDF5");
    GeometryConfiguration *gconf = NULL;
    gconf = GeometryConfiguration::from_legacy(data_file);
    int s_nsph = sconf->number_of_spheres;
    int nsph = gconf->number_of_spheres;
    int configurations = sconf->configurations;
    logger->log(" done.\n");
    // Sanity check on number of sphere consistency, should always be verified
    if (s_nsph == nsph) {
      // Memory requirements of MPI process 0, which also stores the results
      // of all the scales: the full output, one iteration data instance per
      // thread and the one-scale output of each thread other than thread 0.
      int omp_wavelength_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
      {
	if (omp_get_thread_num() == 0) omp_wavelength_threads = omp_get_num_threads();
      }
#endif // _OPENMP
      const long sid_size_bytes = SphereIterationData::get_size(gconf, sconf);
      const long output_size_bytes = SphereOutputInfo::get_size(sconf, gconf);
      const long scale_output_bytes = SphereOutputInfo::get_size(sconf, gconf, 1, 1);
      const long requested_ram_bytes = output_size_bytes
	+ omp_wavelength_threads * sid_size_bytes
	+ (omp_wavelength_threads - 1) * scale_output_bytes;
      char ram_line[128];
      sprintf(
	ram_line, "output data require %.5lg GiB of RAM.\n",
	output_size_bytes / 1024.0 / 1024.0 / 1024.0
      );
      logger->log(ram_line);
      sprintf(
	ram_line, "code execution needs %.5lg GiB of RAM (%d thread(s)).\n",
	requested_ram_bytes / 1024.0 / 1024.0 / 1024.0, omp_wavelength_threads
      );
      logger->log(ram_line);
      // Needs of all the processes, against the memory of their nodes (see
      // resources.h).
      MemoryDemand demand;
      demand.full_output = output_size_bytes;
      demand.scale_output = scale_output_bytes;
      demand.per_thread = sid_size_bytes;
      demand.device_per_thread = 0;
      check_memory(survey, demand, gconf->host_ram_gb, 0.0, gconf->memory_policy, logger);
      ScatteringAngles *p_sa = new ScatteringAngles(gconf);
      SphereIterationData *sid = new SphereIterationData(gconf, sconf, mpidata, 0);
      SphereOutputInfo *p_output = new SphereOutputInfo(sconf, gconf, mpidata);
      // FILE *output = fopen((output_path + "/c_OSPH").c_str(), "w");
      const double half_pi = acos(0.0);
      const double pi = 2.0 * half_pi;
      sid->c1->gcs = 0.0;
      for (int i116 = 0; i116 < nsph; i116++) {
	int i = i116 + 1;
	int iogi = sid->c1->iog[i116];
	if (iogi >= i) {
	  double gcss = pi * sid->c1->ros[i116] * sid->c1->ros[i116];
	  sid->c1->gcsv[i116] = gcss;
	  int nsh = sid->c1->nshl[i116];
	  for (int j115 = 0; j115 < nsh; j115++) {
	    sid->c1->rc[i116][j115] = sconf->get_rcf(i116, j115) * sid->c1->ros[i116];
	  }
	}
	sid->c1->gcs += sid->c1->gcsv[iogi - 1];
      }
      thdps(gconf->l_max, sid->zpv);
      double exdc = sconf->exdc;
      double exri = sqrt(exdc);

      int imode = 10;

      if (sconf->idfc < 0) {
	sid->vk = sid->xip * sid->wn;
	p_output->vec_vk[0] = sid->vk;
      }

      int jer = 0;

      share_setup_status(mpidata, 0);

      setup_shared = true;

      // here go the calls that send data to be duplicated on other MPI processes from process 0 to others, using MPI broadcasts, but only if MPI is actually used
#ifdef MPI_VERSION
      if (mpidata->mpirunning) {
	gconf->mpibcast(mpidata);
	sconf->mpibcast(mpidata);	    
	sid->mpibcast(mpidata);
	p_sa->mpibcast(mpidata);
	int send_barrier = MPI_Barrier(MPI_COMM_WORLD);
      }
#endif
      // Create this variable and initialise it with a default here, so that it is defined anyway, with or without OpenMP support enabled
      int ompnumthreads = 1;
      // this is for MPI process 0 (or even if we are not using MPI at all)
      int myjxi488startoffset = 0;
      int myMPIstride = ompnumthreads;
      int myMPIblock = ompnumthreads;
      SphereOutputInfo **p_outarray = NULL;
      // Errors of the threads in the loop on scales (see ScaleLoopErrors).
      ScaleLoopErrors loop_errors;
      // T-matrix of scale JWTM, if a thread of this process computes it.
      TransitionMatrix *tmatrix = NULL;
      int loop_status = 0;

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
	  p_outarray = new SphereOutputInfo*[ompnumthreads]();
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
	SphereIterationData *sid_2 = NULL;
	SphereOutputInfo *p_output_2 = NULL;
	// for threads other than the 0, create distinct copies of all relevant data, while for thread 0 just define new references / pointers to the original ones
	if (myompthread == 0) {
	  sid_2 = sid;
	  // OMP thread 0 of MPI process 0 holds the pointer to the full output structure
	  p_output_2 = p_output;
	  p_outarray[0] = p_output_2;
	} else {
	  // this is not thread 0, so do create fresh copies of all local variables
	  try {
	    sid_2 = new SphereIterationData(*sid);
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
	for (int ixi488=1; ixi488<=sid->number_of_scales; ixi488 +=myMPIstride) {
	  // the parallel loop over MPI processes covers a different set of indices for each thread
#pragma omp barrier
	  int myjxi488 = ixi488+myompthread;
#pragma omp barrier

	  // each MPI process handles a number of contiguous scales corresponding to its number of OMP threads at this omp level of parallelism
	  // A failing thread records its error, and all the threads stop after
	  // this block of scales (see ScaleLoopErrors).
	  if (sid_2 != NULL) {
	    try {
	      if (myjxi488 <= sid->number_of_scales) {
		if (myompthread > 0) {
		  // UPDATE: non-0 threads need to allocate memory for one scale at a time.
		  p_output_2 = new SphereOutputInfo(sconf, gconf, mpidata, myjxi488, 1);
		  p_outarray[myompthread] = p_output_2;
		}
		int jer = sphere_jxi488_cycle(myjxi488 - 1, sconf, gconf, p_sa, sid_2, p_output_2, output_path, tmatrix);
	      } else {
		if (myompthread > 0) {
		  // If there is no input for this thread, mark to skip.
		  p_outarray[myompthread] = new SphereOutputInfo(1);
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

#ifdef MPI_VERSION
	    if (mpidata->mpirunning) {
	      // only go through this if MPI has been actually used
	      for (int rr=1; rr<mpidata->nprocs; rr++) {
		// get the data from process rr by receiving it in total memory structure
		p_outarray[0]->mpireceive(mpidata, rr);
		
		int test = MPI_Barrier(MPI_COMM_WORLD);
	      }
	    }
#endif
	  }
	  // end block writing to disk
#pragma omp barrier

	} // ixi488 strided MPI loop
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
	delete sid_2;
      } // OMP parallel
      delete p_sa;
      // Stop here if a thread of any process failed: no partial output files.
      loop_errors.rethrow(loop_status);
      // Only MPI process 0 writes the T-matrix (all the processes take part).
      TransitionMatrix::write_from_main_process(mpidata, tmatrix, output_path);
      delete tmatrix;
      check_output_file(output_path + "/c_OSPH.hd5", "output file");
      check_output_file(output_path + "/c_OSPH", "output file");
      p_output->write(output_path + "/c_OSPH.hd5", "HDF5");
      p_output->write(output_path + "/c_OSPH", "LEGACY");
      delete p_output;
      logger->log("Finished. Output written to " + output_path + "/c_OSPH.\n");
    } else { // NSPH mismatch between geometry and scatterer configurations.
      throw UnrecognizedConfigurationException(
        "inconsistent configurations: the scatterer configuration describes "
	+ to_string(s_nsph) + " spheres, the geometry configuration " + to_string(nsph)
      );
    }
    delete sconf;
    delete gconf;
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
  else { // Instruction block for MPI processes other than 0.
    // here go the code for MPI processes other than 0
    // copy gconf, sconf, cid and p_scattering_angles from MPI process 0
    receive_setup_status(mpidata);
    GeometryConfiguration *gconf = new GeometryConfiguration(mpidata);
    ScattererConfiguration *sconf = new ScattererConfiguration(mpidata);
    SphereIterationData *sid = new SphereIterationData(mpidata, device_count);
    ScatteringAngles *p_sa = new ScatteringAngles(mpidata);
    int recv_barrier = MPI_Barrier(MPI_COMM_WORLD);
    
    // Create this variable and initialise it with a default here, so that it is defined anyway, with or without OpenMP support enabled
    int ompnumthreads = 1;
    SphereOutputInfo **p_outarray = NULL;
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
	p_outarray = new SphereOutputInfo*[ompnumthreads]();
      }
#pragma omp barrier
      // To test parallelism, I will now start feeding this function with "clean" copies of the parameters, so that they will not be changed by previous iterations, and each one will behave as the first one. Define all (empty) variables here, so they have the correct scope, then they get different definitions depending on thread number
      SphereIterationData *sid_2 = NULL;
      SphereOutputInfo *p_output_2 = NULL;
      // PLACEHOLDER
      // for threads other than the 0, create distinct copies of all relevant data, while for thread 0 just define new references / pointers to the original ones
      if (myompthread == 0) {
	sid_2 = sid;
      } else {
	// this is not thread 0, so do create fresh copies of all local variables
	try {
	  sid_2 = new SphereIterationData(*sid);
	} catch (...) {
	  loop_errors.capture();
	}
      }
      // make sure all threads align here: I don't want the following loop to accidentally start for thread 0, possibly modifying some variables before they are copied by all other threads
#pragma omp barrier
      // ok, now I can actually start the parallel calculations
      for (int ixi488=1; ixi488<=sid->number_of_scales; ixi488 +=myMPIstride) {
	// the parallel loop over MPI processes covers a different set of indices for each thread
#pragma omp barrier
	int myjxi488 = ixi488 + myjxi488startoffset + myompthread;
	// each thread opens new virtual files and stores their pointers in the shared array
#pragma omp barrier
	if (myompthread==0) logger->log("Syncing OpenMP threads and starting the loop on wavelengths\n");
	// ok, now I can actually start the parallel calculations
	// each MPI process handles a number of contiguous scales corresponding to its number of OMP threads at this omp level of parallelism
	// A failing thread records its error, and all the threads stop after
	// this block of scales (see ScaleLoopErrors).
	if (sid_2 != NULL) {
	  try {
	    if (myjxi488 <= sid->number_of_scales) {
	      if (myompthread > 0) {
		// UPDATE: non-0 threads need to allocate memory for one scale at a time.
		p_output_2 = new SphereOutputInfo(sconf, gconf, mpidata, myjxi488, 1);
		p_outarray[myompthread] = p_output_2;
	      } else {
		// Thread 0 of non-zero MPI processes needs to allocate memory for the
		// output of all threads _doing something_.
		int iterstodo = sid->number_of_scales - myjxi488 + 1;
		if (iterstodo > ompnumthreads) iterstodo = ompnumthreads;
		p_output_2 = new SphereOutputInfo(sconf, gconf, mpidata, myjxi488, iterstodo);
		p_outarray[0] = p_output_2;
	      }
	      int jer = sphere_jxi488_cycle(myjxi488 - 1, sconf, gconf, p_sa, sid_2, p_output_2, output_path, tmatrix);
	    } else {
	      p_outarray[myompthread] = new SphereOutputInfo(1);
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
      delete sid_2;
    } // OMP parallel
    delete p_sa;
    delete sconf;
    delete gconf;
    // Stop here if a thread of any process failed.
    loop_errors.rethrow(loop_status);
    // Only MPI process 0 writes the T-matrix (all the processes take part).
    TransitionMatrix::write_from_main_process(mpidata, tmatrix, output_path);
    delete tmatrix;
  } // End instructions block for MPI non-0 processes.
#endif // MPI_VERSION
  delete logger;
} // sphere()

int sphere_jxi488_cycle(
  int jxi488, ScattererConfiguration *sconf, GeometryConfiguration *gconf,
  ScatteringAngles *sa, SphereIterationData *sid, SphereOutputInfo *oi,
  const string& output_path, TransitionMatrix *&tmatrix
) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  const double half_pi = acos(0.0);
  const double pi = 2.0 * half_pi;
  int jer = 0;
  Logger *logger = new Logger(LOG_INFO);
  string message = "INIT";
  int oindex = 0;
  int jxi = jxi488 + 1;
  int jxindex = jxi - oi->first_xi;
  // this is now part of sid, so don't mess with it here, just copy it by reference
  bool &is_first_scale = sid->is_first_scale;
  //bool is_first_scale = (jxi == 1);
  int nsph = gconf->number_of_spheres;
  //int l_max = gconf->l_max;
  int in_pol = gconf->in_pol;
  int npnt = gconf->npnt;
  int npntts = gconf->npntts;
  int jwtm = gconf->jwtm;
  double wn = sconf->wp / 3.0e8;
  double sqsfi = 1.0;
  int idfc = sconf->idfc;
  int nxi = sconf->number_of_scales;
  int configurations = sconf->configurations;
  int ndirs = sa->nkks;
  int lcalc;
  int jw, isq, ibf;
  logger->log("running scale iteration " + to_string(jxi) + " of " + to_string(nxi) + ".\n");
  double vk, vkarg;
  double xi = sconf->get_scale(jxi488);
  double exdc = sconf->exdc;
  double exri = sqrt(exdc);
  if (idfc >= 0) {
    vk = xi * wn;
    vkarg = vk;
    oi->vec_vk[jxindex] = vk;
    oi->vec_xi[jxindex] = xi;
  } else { // IDFC < 0
    vk = sconf->xip * wn;
    vkarg = xi * vk;
    sqsfi = 1.0 / (xi * xi);
    oi->vec_vk[jxindex] = vk;
    oi->vec_xi[jxindex] = xi;
  }
  // Dynamic order check
  const int max_lm = gconf->l_max;
  int l_max = gconf->l_max;
  double xi_factor = (sconf->reference_variable_name.compare("XIV") == 0) ? xi : 1.0;
  const double alamb = 2.0 * pi * xi_factor / vk;
  double size_par_lm = 2.0 * pi * sqrt(exdc) * sconf->get_max_radius() / alamb;
  int recommended_lm = 4 + (int)ceil(size_par_lm + 4.05 * pow(size_par_lm, 1.0 / 3.0));
  if (recommended_lm != l_max) {
    if (recommended_lm < max_lm) {
      if (gconf->dyn_order_flag > 0) {
	int new_lm = recommended_lm;
	message = "lowering internal order from " + to_string(max_lm) + " to "
	  + to_string(recommended_lm) + " for scale iteration " + to_string(jxi488 + 1) + ".\n";
	logger->log(message, LOG_INFO);
	sid->update_order(new_lm);
	is_first_scale = true;
	// jw = 1;
	l_max = new_lm;
      } else {
	message = "internal order " + to_string(max_lm) + " for scale iteration "
	  + to_string(jxi488 + 1) + " too high (recommended order is " + to_string(recommended_lm)
	  + ").\n";
      logger->log(message, LOG_WARN);
      }
    }
  }
  // End of dynamic order check
  double thsca = (gconf->isam > 1) ? sa->ths - sa->th : 0.0;
  for (int i132 = 0; i132 < nsph; i132++) {
    int i = i132 + 1;
    int iogi = sid->c1->iog[i132];
    if (iogi != i) {
      for (int l123 = 0; l123 < l_max; l123++) {
	sid->c1->rmi[l123][i132] = sid->c1->rmi[l123][iogi - 1];
	sid->c1->rei[l123][i132] = sid->c1->rei[l123][iogi - 1];
      }
      continue; // i132
    }
    // label 125
    int nsh = sid->c1->nshl[i132];
    int ici = (nsh + 1) / 2;
    if (idfc == 0) {
      for (int ic = 0; ic < ici; ic++)
	sid->c1->dc0[ic] = sconf->get_dielectric_constant(ic, i132, jxi488); // WARNING: IDFC=0 is not tested!
    } else { // IDFC != 0
      if (is_first_scale) {
	int scale_to_extract = (sconf->reference_variable_name.compare("XIV") == 0) ? 0 : jxi488;
	for (int ic = 0; ic < ici; ic++) {
	  sid->c1->dc0[ic] = sconf->get_dielectric_constant(ic, i132, scale_to_extract);
	}
      }
    }
    if (nsh % 2 == 0) sid->c1->dc0[ici] = exdc;
    dme(l_max, i, npnt, npntts, vkarg, exdc, exri, sid->c1, jer, lcalc, sid->arg);
    if (jer != 0) {
      oi->vec_ier[jxindex] = 1;
      oi->lcalc = lcalc;
      delete logger;
      return jer;
    }
  } // i132
  if (idfc >= 0 and nsph == 1 and jxi == jwtm) {
    // This is the condition that produces the transition matrix: keep it,
    // MPI process 0 writes it after the loop on scales.
    tmatrix = new TransitionMatrix(
      l_max, vk, exri, sid->c1->rmi, sid->c1->rei, sconf->get_type_radius(0)
    );
  }
  double cs0 = 0.25 * vk * vk * vk / half_pi;
  sscr0(sid->tfsas, nsph, l_max, vk, exri, sid->c1);
  double sqk = vk * vk * exdc;
  aps(sid->zpv, l_max, nsph, sid->c1, sqk, sid->gaps);
  rabas(in_pol, l_max, nsph, sid->c1, sid->tqse, sid->tqspe, sid->tqss, sid->tqsps);
  int last_configuration = 0;
  for (int i170 = 0; i170 < nsph; i170++) {
    int i = i170 + 1;
    if (sid->c1->iog[i170] >= i) {
      last_configuration++;
      oindex = jxindex * configurations + last_configuration - 1;
      double albeds = sid->c1->sscs[i170] / sid->c1->sexs[i170];
      sid->c1->sqscs[i170] *= sqsfi;
      sid->c1->sqabs[i170] *= sqsfi;
      sid->c1->sqexs[i170] *= sqsfi;
      if (sid->c1->nshl[i170] != 1) {
	oi->vec_sphere_ref_indices[oindex] = cc0;
	oi->vec_sphere_sizes[oindex] = sid->c1->vsz[i170];
      } else {
	oi->vec_sphere_ref_indices[oindex] = sid->c1->vkt[i170];
	oi->vec_sphere_sizes[oindex] = sid->c1->vsz[i170];
      }
      oi->vec_scs[oindex] = sid->c1->sscs[i170];
      oi->vec_abs[oindex] = sid->c1->sabs[i170];
      oi->vec_exs[oindex] = sid->c1->sexs[i170];
      oi->vec_albeds[oindex] = albeds;
      oi->vec_scsrt[oindex] = sid->c1->sqscs[i170];
      oi->vec_absrt[oindex] = sid->c1->sqabs[i170];
      oi->vec_exsrt[oindex] = sid->c1->sqexs[i170];
      oi->vec_fsas[oindex] = sid->c1->fsas[i170];
      double csch = 2.0 * vk * sqsfi / sid->c1->gcsv[i170];
      dcomplex s0 = sid->c1->fsas[i170] * exri;
      double qschu = csch * imag(s0);
      double pschu = csch * real(s0);
      double s0mag = cs0 * cabs(s0);
      oi->vec_qschu[oindex] = qschu;
      oi->vec_pschu[oindex] = pschu;
      oi->vec_s0mag[oindex] = s0mag;
      double rapr = sid->c1->sexs[i170] - sid->gaps[i170];
      double cosav = sid->gaps[i170] / sid->c1->sscs[i170];
      oi->vec_cosav[oindex] = cosav;
      oi->vec_raprs[oindex] = rapr;
      oi->vec_tqek1[oindex] = sid->tqse[0][i170];
      oi->vec_tqsk1[oindex] = sid->tqss[0][i170];
      oi->vec_tqek2[oindex] = sid->tqse[1][i170];
      oi->vec_tqsk2[oindex] = sid->tqss[1][i170];
    } // End if iog[i170] >= i
  } // i170 loop
  if (nsph != 1) {
    oi->vec_fsat[jxindex] = sid->tfsas;
    double csch = 2.0 * vk * sqsfi / sid->c1->gcs;
    dcomplex s0 = sid->tfsas * exri;
    double qschu = csch * imag(s0);
    double pschu = csch * real(s0);
    double s0mag = cs0 * cabs(s0);
    oi->vec_qschut[jxindex] = qschu;
    oi->vec_pschut[jxindex] = pschu;
    oi->vec_s0magt[jxindex] = s0mag;
  }
  double th = sa->th;
  int done_dirs = 0;
  int nks = sa->nths * sa->nphs;
  int nkks = sa->nth * sa->nph * nks;
  int nth = sa->nth;
  int nph = sa->nph;
  // With ISAM < 0 the forces are not computed, but printed: as in the
  // FORTRAN (whose variables start at zero), they are 0.
  double frx = 0.0, fry = 0.0, frz = 0.0;
  for (int jth486 = 0; jth486 < sa->nth; jth486++) { // OpenMP parallelizable section
    int jth = jth486 + 1;
    double ph = sa->ph;
    for (int jph484 = 0; jph484 < sa->nph; jph484++) {
      int jph = jph484 + 1;
      bool goto182 = (sa->nk == 1) && (!is_first_scale);
      if (!goto182) {
	upvmp(th, ph, 0, sid->cost, sid->sint, sid->cosp, sid->sinp, sid->u, sid->upmp, sid->unmp);
      }
      if (gconf->isam >= 0) {
	wmamp(0, sid->cost, sid->sint, sid->cosp, sid->sinp, in_pol, l_max, 0, nsph, sid->argi, sid->u, sid->upmp, sid->unmp, sid->c1);
	for (int i183 = 0; i183 < nsph; i183++) {
	  double rapr = sid->c1->sexs[i183] - sid->gaps[i183];
	  frx = rapr * sid->u[0];
	  fry = rapr * sid->u[1];
	  frz = rapr * sid->u[2];
	}
	jw = 1;
      }
      double thsl = sa->ths;
      double phsph = 0.0;
      for (int jths482 = 0; jths482 < sa->nths; jths482++) {
	int jths = jths482 + 1;
	double ths = thsl;
	int icspnv = 0;
	if (gconf->isam > 1) ths = th + thsca;
	if (gconf->isam >= 1) {
	  phsph = 0.0;
	  if ((ths < 0.0) || (ths > 180.0)) phsph = 180.0;
	  if (ths < 0.0) ths *= -1.0;
	  if (ths > 180.0) ths = 360.0 - ths;
	  if (phsph != 0.0) icspnv = 1;
	}
	double phs = sa->phs;
	for (int jphs480 = 0; jphs480 < sa->nphs; jphs480++) {
	  int jphs = jphs480 + 1;
	  if (gconf->isam >= 1) {
	    phs = ph + phsph;
	    if (phs >= 360.0) phs -= 360.0;
	  }
	  // label 188: the scattered direction is computed again unless it
	  // cannot change, i.e. it is fixed (ISAM < 1) or there is only one
	  // incident direction (with ISAM >= 1 it follows the incident one).
	  // This differs from the original FORTRAN code and is under evaluation
	  // by its authors: see dev_notes/questions-for-authors.md, question 1.
	  bool goto190 = (nks == 1) && (gconf->isam < 1 || sa->nth * sa->nph == 1)
	    && ((!is_first_scale) || (jth > 1) || (jph > 1));
	  if (!goto190) {
	    upvmp(ths, phs, icspnv, sid->costs, sid->sints, sid->cosps, sid->sinps, sid->us, sid->upsmp, sid->unsmp);
	    if (gconf->isam >= 0) {
	      wmamp(2, sid->costs, sid->sints, sid->cosps, sid->sinps, in_pol, l_max, 0, nsph, sid->args, sid->us, sid->upsmp, sid->unsmp, sid->c1);
	    }
	  }
	  // label 190 (FORTRAN: skipped if NKKS = 1 after the first scale)
	  if (nkks != 1 || is_first_scale) {
	    upvsp(
	      sid->u, sid->upmp, sid->unmp, sid->us, sid->upsmp, sid->unsmp,
	      sid->up, sid->un, sid->ups, sid->uns, sid->duk, isq, ibf,
	      sid->scan, sid->cfmp, sid->sfmp, sid->cfsp, sid->sfsp
	    );
	    if (gconf->isam < 0) {
	      // the scattering plane vectors of upvsp() are kept (label 194)
	      wmasp(
		sid->cost, sid->sint, sid->cosp, sid->sinp, sid->costs,
		sid->sints, sid->cosps, sid->sinps, sid->u, sid->up,
		sid->un, sid->us, sid->ups, sid->uns, isq, ibf, in_pol,
		l_max, 0, nsph, sid->argi, sid->args, sid->c1
	      );
	    } else { // label 192
	      for (int i193 = 0; i193 < 3; i193++) {
		sid->un[i193] = sid->unmp[i193];
		sid->uns[i193] = sid->unsmp[i193];
	      }
	    }
	  }
	  if (gconf->isam < 0) jw = 1;
	  if (jw != 0) {
	    jw = 0;
	  }
	  oi->vec_dir_scand[done_dirs] = sid->scan;
	  oi->vec_dir_cfmp[done_dirs] = sid->cfmp;
	  oi->vec_dir_cfsp[done_dirs] = sid->cfsp;
	  oi->vec_dir_sfmp[done_dirs] = sid->sfmp;
	  oi->vec_dir_sfsp[done_dirs] = sid->sfsp;
	  if (gconf->isam >= 0) {
	    oi->vec_dir_un[3 * done_dirs] = sid->un[0];
	    oi->vec_dir_un[3 * done_dirs + 1] = sid->un[1];
	    oi->vec_dir_un[3 * done_dirs + 2] = sid->un[2];
	    oi->vec_dir_uns[3 * done_dirs] = sid->uns[0];
	    oi->vec_dir_uns[3 * done_dirs + 1] = sid->uns[1];
	    oi->vec_dir_uns[3 * done_dirs + 2] = sid->uns[2];
	  } else {
	    oi->vec_dir_un[3 * done_dirs] = sid->un[0];
	    oi->vec_dir_un[3 * done_dirs + 1] = sid->un[1];
	    oi->vec_dir_un[3 * done_dirs + 2] = sid->un[2];
	    oi->vec_dir_uns[3 * done_dirs] = 0.0;
	    oi->vec_dir_uns[3 * done_dirs + 1] = 0.0;
	    oi->vec_dir_uns[3 * done_dirs + 2] = 0.0;
	  }
	  sscr2(nsph, l_max, vk, exri, sid->c1);
	  last_configuration = 0;
	  for (int ns226 = 0; ns226 < nsph; ns226++) {
	    int ns = ns226 + 1;
	    oindex = jxindex * nsph * ndirs + nsph * done_dirs + ns226;
	    oi->vec_dir_sas11[oindex] = sid->c1->sas[ns226][0][0];
	    oi->vec_dir_sas21[oindex] = sid->c1->sas[ns226][1][0];
	    oi->vec_dir_sas12[oindex] = sid->c1->sas[ns226][0][1];
	    oi->vec_dir_sas22[oindex] = sid->c1->sas[ns226][1][1];
	    oi->vec_dir_fx[jxindex * nsph * nth * nph + nsph * nph * (jth - 1) + nsph * (jph - 1) + ns226] = frx;
	    oi->vec_dir_fy[jxindex * nsph * nth * nph + nsph * nph * (jth - 1) + nsph * (jph - 1) + ns226] = fry;
	    oi->vec_dir_fz[jxindex * nsph * nth * nph + nsph * nph * (jth - 1) + nsph * (jph - 1) + ns226] = frz;
	    for (int i225 = 0; i225 < 16; i225++) sid->c1->vint[i225] = sid->c1->vints[ns226][i225];
	    mmulc(sid->c1->vint, sid->cmullr, sid->cmul);
	    for (int imul = 0; imul < 4; imul++) {
	      int muls_index = 16 * jxindex * nsph * ndirs + 16 * nsph * done_dirs + 4 * imul;
	      for (int jmul = 0; jmul < 4; jmul++) {
		oi->vec_dir_muls[muls_index + jmul] = sid->cmul[imul][jmul];
	      }
	    }
	    for (int imul = 0; imul < 4; imul++) {
	      int muls_index = 16 * jxindex * nsph * ndirs + 16 * nsph * done_dirs + 4 * imul;
	      for (int jmul = 0; jmul < 4; jmul++) {
		oi->vec_dir_mulslr[muls_index + jmul] = sid->cmullr[imul][jmul];
	      }
	    }
	  } // ns226 loop
	  if (gconf->isam < 1) phs += sa->phsstp;
	  done_dirs++;
	} // jphs480 loop
	if (gconf->isam <= 1) thsl += sa->thsstp;
      } // jths482 loop
      ph += sa->phstp;
    } // jph484 loop on elevation
    th += sa->thstp;
  } // jth486 loop on azimuth
  // at the end of this iteration make sure to set is_first_scale to false, the next iteration will reinitialise _only_ if the order changes
  is_first_scale = 0;
  oi->vec_jxi[jxindex] = jxi;
  logger->log("finished scale iteration " + to_string(jxi) + " of " + to_string(nxi) + ".\n");
  delete logger;
  return jer;
}

// >>> IMPLEMENTATION OF SphereIterationData CLASS <<<
SphereIterationData::SphereIterationData(
  const GeometryConfiguration *gconf, const ScattererConfiguration *sconf,
  const mixMPI *mpidata, const int device_count
) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  _nsph = gconf->number_of_spheres;
  _lm = gconf->l_max;
  arg = cc0;
  s0 = cc0;
  tfsas = cc0;
  c1 = new ParticleDescriptorSphere(gconf, sconf);
  argi = new double[1]();
  args = new double[1]();
  cost = 0.0;
  sint = 0.0;
  cosp = 0.0;
  sinp = 0.0;
  costs = 0.0;
  sints = 0.0;
  cosps = 0.0;
  sinps = 0.0;
  scan = 0.0;
  cfmp = 0.0;
  cfsp = 0.0;
  sfmp = 0.0;
  sfsp = 0.0;
  wn = sconf->wp / 3.0e8;
  xip = sconf->xip;
  vk = 0.0;
  // Scale block initialization
  number_of_scales = sconf->number_of_scales;
  xiblock = (int) ceil(((double) (sconf->number_of_scales-1))/((double) mpidata->nprocs));
  lastxi = ((mpidata->rank+1) * xiblock)+1;
  firstxi = lastxi-xiblock+1;
  if (lastxi > sconf->number_of_scales) lastxi = sconf->number_of_scales;
  // End of scale block initialization
  gaps = new double[_nsph]();
  duk = new double[3]();
  u = new double[3]();
  us = new double[3]();
  un = new double[3]();
  uns = new double[3]();
  up = new double[3]();
  ups = new double[3]();
  upmp = new double[3]();
  upsmp = new double[3]();
  unmp = new double[3]();
  unsmp = new double[3]();
  vec_cmul = new double[16]();
  vec_cmullr = new double[16]();
  cmul = new double*[4];
  cmullr = new double*[4];
  for (int ci = 0; ci < 4; ci++) {
    cmul[ci] = (vec_cmul + 4 * ci);
    cmullr[ci] = (vec_cmullr + 4 * ci);
  }
  vec_tqspe = new dcomplex[2 * _nsph]();
  vec_tqsps = new dcomplex[2 * _nsph]();
  vec_tqse = new double[2 * _nsph]();
  vec_tqss = new double[2 * _nsph]();
  tqspe = new dcomplex*[2];
  tqsps = new dcomplex*[2];
  tqse = new double*[2];
  tqss = new double*[2];
  for (int ti = 0; ti < 2; ti++) {
    tqspe[ti] = (vec_tqspe + _nsph * ti);
    tqsps[ti] = (vec_tqsps + _nsph * ti);
    tqse[ti] = (vec_tqse + _nsph * ti);
    tqss[ti] = (vec_tqss + _nsph * ti);
  }
  vec_zpv = new double[_lm * 12]();
  zpv = new double***[_lm];
  for (int zi = 0; zi < _lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      int vec_index = 12 * zi + 4 * zj;
      zpv[zi][zj] = new double*[2];
      zpv[zi][zj][0] = (vec_zpv + vec_index);
      zpv[zi][zj][1] = (vec_zpv + vec_index + 2);
    }
  }
  is_first_scale = 1;
}

long SphereIterationData::get_size(const GeometryConfiguration *gconf, const ScattererConfiguration *sconf) {
  const long nsph = gconf->number_of_spheres;
  const long lm = gconf->l_max;
  long result = sizeof(SphereIterationData);
  result += ParticleDescriptorSphere::get_size(gconf, sconf); // c1
  result += 2 * sizeof(double); // argi, args
  result += sizeof(double) * nsph; // gaps
  result += 11 * sizeof(double) * 3; // duk, u, us, un, uns, up, ups, upmp, upsmp, unmp, unsmp
  result += 2 * sizeof(double) * 16; // vec_cmul, vec_cmullr
  result += 2 * sizeof(double *) * 4; // cmul, cmullr
  result += 2 * sizeof(dcomplex) * 2 * nsph; // vec_tqspe, vec_tqsps
  result += 2 * sizeof(double) * 2 * nsph; // vec_tqse, vec_tqss
  result += 2 * sizeof(dcomplex *) * 2 + 2 * sizeof(double *) * 2; // tqspe, tqsps, tqse, tqss
  result += sizeof(double) * lm * 12; // vec_zpv
  result += sizeof(double ***) * lm; // zpv
  result += sizeof(double **) * 3 * lm; // zpv[i]
  result += sizeof(double *) * 2 * 3 * lm; // zpv[i][j]
  return result;
}

SphereIterationData::SphereIterationData(const SphereIterationData &rhs) {
  _nsph = rhs._nsph;
  _lm = rhs._lm;
  arg = rhs.arg;
  s0 = rhs.s0;
  tfsas = rhs.tfsas;
  c1 = new ParticleDescriptorSphere(reinterpret_cast<ParticleDescriptorSphere &>(*(rhs.c1)));
  argi = new double[1];
  args = new double[1];
  argi[0] = rhs.argi[0];
  args[0] = rhs.args[0];
  cost = rhs.cost;
  sint = rhs.sint;
  cosp = rhs.cosp;
  sinp = rhs.sinp;
  costs = rhs.costs;
  sints = rhs.sints;
  cosps = rhs.cosps;
  sinps = rhs.sinps;
  scan = rhs.scan;
  cfmp = rhs.cfmp;
  cfsp = rhs.cfsp;
  sfmp = rhs.sfmp;
  sfsp = rhs.sfsp;
  wn = rhs.wn;
  xip = rhs.xip;
  vk = rhs.vk;
  // Scale block initialization
  number_of_scales = rhs.number_of_scales;
  xiblock = rhs.xiblock;
  lastxi = rhs.lastxi;
  firstxi = rhs.firstxi;
  // End of scale block initialization
  gaps = new double[_nsph];
  for (int si = 0; si < _nsph; si++) gaps[si] = rhs.gaps[si];
  duk = new double[3];
  u = new double[3];
  us = new double[3];
  un = new double[3];
  uns = new double[3];
  up = new double[3];
  ups = new double[3];
  upmp = new double[3];
  upsmp = new double[3];
  unmp = new double[3];
  unsmp = new double[3];
  for (int ui = 0; ui < 3; ui++) {
    duk[ui] = rhs.duk[ui];
    u[ui] = rhs.u[ui];
    us[ui] = rhs.us[ui];
    un[ui] = rhs.un[ui];
    uns[ui] = rhs.us[ui];
    up[ui] = rhs.up[ui];
    ups[ui] = rhs.ups[ui];
    upmp[ui] = rhs.upmp[ui];
    upsmp[ui] = rhs.upsmp[ui];
    unmp[ui] = rhs.unmp[ui];
    unsmp[ui] = rhs.unsmp[ui];
  }
  vec_cmul = new double[16];
  vec_cmullr = new double[16];
  for (int mi = 0; mi < 16; mi++) {
    vec_cmul[mi] = rhs.vec_cmul[mi];
    vec_cmullr[mi] = rhs.vec_cmullr[mi];
  }
  cmul = new double*[4];
  cmullr = new double*[4];
  for (int ci = 0; ci < 4; ci++) {
    cmul[ci] = (vec_cmul + 4 * ci);
    cmullr[ci] = (vec_cmullr + 4 * ci);
  }
  vec_tqspe = new dcomplex[2 * _nsph]();
  vec_tqsps = new dcomplex[2 * _nsph]();
  vec_tqse = new double[2 * _nsph]();
  vec_tqss = new double[2 * _nsph]();
  for (int ti = 0; ti < 2 * _nsph; ti++) {
    vec_tqspe[ti] = rhs.vec_tqspe[ti];
    vec_tqsps[ti] = rhs.vec_tqsps[ti];
    vec_tqse[ti] = rhs.vec_tqse[ti];
    vec_tqss[ti] = rhs.vec_tqss[ti];
  }
  tqspe = new dcomplex*[2];
  tqsps = new dcomplex*[2];
  tqse = new double*[2];
  tqss = new double*[2];
  for (int ti = 0; ti < 2; ti++) {
    tqspe[ti] = (vec_tqspe + _nsph * ti);
    tqsps[ti] = (vec_tqsps + _nsph * ti);
    tqse[ti] = (vec_tqse + _nsph * ti);
    tqss[ti] = (vec_tqss + _nsph * ti);
  }
  vec_zpv = new double[_lm * 12];
  for (int zi = 0; zi < _lm * 12; zi++) vec_zpv[zi] = rhs.vec_zpv[zi];
  zpv = new double***[_lm];
  for (int zi = 0; zi < _lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      int vec_index = 12 * zi + 4 * zj;
      zpv[zi][zj] = new double*[2];
      zpv[zi][zj][0] = (vec_zpv + vec_index);
      zpv[zi][zj][1] = (vec_zpv + vec_index + 2);
    }
  }
  is_first_scale = rhs.is_first_scale;
}

SphereIterationData::~SphereIterationData() {
  int lm = c1->li;
  delete c1;
  delete[] argi;
  delete[] args;
  delete[] gaps;
  delete[] duk;
  delete[] u;
  delete[] us;
  delete[] un;
  delete[] uns;
  delete[] up;
  delete[] ups;
  delete[] upmp;
  delete[] upsmp;
  delete[] unmp;
  delete[] unsmp;
  delete[] vec_cmul;
  delete[] vec_cmullr;
  delete[] cmul;
  delete[] cmullr;
  delete[] vec_tqspe;
  delete[] vec_tqsps;
  delete[] vec_tqse;
  delete[] vec_tqss;
  delete[] tqspe;
  delete[] tqsps;
  delete[] tqse;
  delete[] tqss;
  delete[] vec_zpv;
  for (int zi = 0; zi < lm; zi++) {
    for (int zj = 0; zj < 3; zj++) {
      delete[] zpv[zi][zj];
    }
    delete[] zpv[zi];
  }
  delete[] zpv;
}

#ifdef MPI_VERSION
SphereIterationData::SphereIterationData(const mixMPI *mpidata, const int device_count) {
  // Buld child object
  c1 = new ParticleDescriptorSphere(mpidata);

  // Collect the scalar values
  MPI_Bcast(&_nsph, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lm, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&arg, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&s0, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tfsas, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sint, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cosp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sinp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&costs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sints, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cosps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sinps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&scan, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&wn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xip, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vk, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xiblock, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&number_of_scales, 1, MPI_INT, 0, MPI_COMM_WORLD);
  lastxi = ((mpidata->rank+1) * xiblock)+1;
  firstxi = lastxi-xiblock+1;
  if (lastxi > number_of_scales) lastxi = number_of_scales;

  // Collect length-1 vectors
  argi = new double[1];
  args = new double[1];
  MPI_Bcast(argi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(args, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Collect vectors whose size depend on NSPH
  gaps = new double[_nsph];
  vec_tqspe = new dcomplex[2 * _nsph];
  vec_tqsps = new dcomplex[2 * _nsph];
  vec_tqse = new double[2 * _nsph];
  vec_tqss = new double[2 * _nsph];
  MPI_Bcast(gaps, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqspe, 2 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqsps, 2 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqse, 2 * _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqss, 2 * _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Collect length-3 vectors
  duk = new double[3];
  u = new double[3];
  us = new double[3];
  un = new double[3];
  uns = new double[3];
  up = new double[3];
  ups = new double[3];
  upmp = new double[3];
  upsmp = new double[3];
  unmp = new double[3];
  unsmp = new double[3];
  MPI_Bcast(duk, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(u, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(us, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(un, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(uns, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(up, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ups, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Collect length-16 vectors
  vec_cmul = new double[16];
  vec_cmullr = new double[16];
  MPI_Bcast(vec_cmul, 16, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_cmullr, 16, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Collect vectors whose size depend on LM
  vec_zpv = new double[12 * _lm];
  MPI_Bcast(vec_zpv, 12 * _lm, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Build the multi-dimensional views on the received vectors, as the
  // regular and copy constructors do.
  cmul = new double*[4];
  cmullr = new double*[4];
  for (int ci = 0; ci < 4; ci++) {
    cmul[ci] = (vec_cmul + 4 * ci);
    cmullr[ci] = (vec_cmullr + 4 * ci);
  }
  tqspe = new dcomplex*[2];
  tqsps = new dcomplex*[2];
  tqse = new double*[2];
  tqss = new double*[2];
  for (int ti = 0; ti < 2; ti++) {
    tqspe[ti] = (vec_tqspe + _nsph * ti);
    tqsps[ti] = (vec_tqsps + _nsph * ti);
    tqse[ti] = (vec_tqse + _nsph * ti);
    tqss[ti] = (vec_tqss + _nsph * ti);
  }
  zpv = new double***[_lm];
  for (int zi = 0; zi < _lm; zi++) {
    zpv[zi] = new double**[3];
    for (int zj = 0; zj < 3; zj++) {
      int vec_index = 12 * zi + 4 * zj;
      zpv[zi][zj] = new double*[2];
      zpv[zi][zj][0] = (vec_zpv + vec_index);
      zpv[zi][zj][1] = (vec_zpv + vec_index + 2);
    }
  }

  MPI_Bcast(&is_first_scale, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

}

int SphereIterationData::mpibcast(const mixMPI *mpidata) {
  int result = 0;
  // Broadcast child object
  c1->mpibcast(mpidata);
  
  // Broadcast scalar values
  MPI_Bcast(&_nsph, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lm, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&arg, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&s0, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tfsas, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sint, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cosp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sinp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&costs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sints, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cosps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sinps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&scan, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sfsp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&wn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xip, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&vk, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xiblock, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&number_of_scales, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Broadcast length-1 vectors
  MPI_Bcast(argi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(args, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Broadcast vectors whose size depend on NSPH
  MPI_Bcast(gaps, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqspe, 2 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqsps, 2 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqse, 2 * _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_tqss, 2 * _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Broadcast length-3 vectors
  MPI_Bcast(duk, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(u, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(us, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(un, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(uns, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(up, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ups, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(upsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(unsmp, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Broadcast length-16 vectors
  MPI_Bcast(vec_cmul, 16, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_cmullr, 16, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Broadcast vectors whose size depend on LM
  MPI_Bcast(vec_zpv, 12 * _lm, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(&is_first_scale, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

  return 0;
}
#endif // MPI_VERSION

int SphereIterationData::update_order(int order) {
  int result = 0;
  int old_lm = _lm;
  if (order != _lm) {
    _lm = order;
    ((ParticleDescriptorSphere *)c1)->update_order(order);
    for (int zi = 0; zi < old_lm; zi++) {
      for (int zj = 0; zj < 3; zj++) {
	delete[] zpv[zi][zj];
      }
      delete[] zpv[zi];
    }
    delete[] zpv;
    delete[] vec_zpv;
    vec_zpv = new double[_lm * 12]();
    zpv = new double***[_lm];
    for (int zi = 0; zi < _lm; zi++) {
      zpv[zi] = new double**[3];
      for (int zj = 0; zj < 3; zj++) {
	int vec_index = 12 * zi + 4 * zj;
	zpv[zi][zj] = new double*[2];
	zpv[zi][zj][0] = (vec_zpv + vec_index);
	zpv[zi][zj][1] = (vec_zpv + vec_index + 2);
      }
    }
    thdps(_lm, zpv);
  }
  return result;
}
// >>> END OF SphereIterationData CLASS IMPLEMENTATION <<<
