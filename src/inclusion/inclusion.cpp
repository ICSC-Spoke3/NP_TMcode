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

/*! \file inclusion.cpp
 *
 * \brief Implementation of the calculation for a sphere with inclusions.
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
#include <cuda_runtime.h>
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

#ifndef INCLUDE_INCLU_SUBS_H_
#include "../include/inclu_subs.h"
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

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

#ifndef INCLUDE_ITERATION_DATA_H_
#include "../include/IterationData.h"
#endif

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
 *  \param cid: `InclusionIterationData *` Pointer to an `InclusionIterationData` object.
 *  \param output: `InclusionOutputInfo *` Pointer to an `InclusionOutputInfo` object.
 *  \param output_path: `const string &` Path to the output directory.
 *  \param vtppoanp: `VirtualBinaryFile *` Pointer to a `VirtualBinaryFile` object.
 */
int inclusion_jxi488_cycle(int jxi488, ScattererConfiguration *sconf, GeometryConfiguration *gconf, ScatteringAngles *sa, InclusionIterationData *cid, InclusionOutputInfo *output, const string& output_path, VirtualBinaryFile *vtppoanp);

/*! \brief C++ implementation of INCLU
 *
 * \param config_file: `string` Name of the configuration file.
 * \param data_file: `string` Name of the input data file.
 * \param output_path: `string` Directory to write the output files in.
 * \param mpidata: `mixMPI *` Pointer to an instance of MPI data settings.
 */
void inclusion(const string& config_file, const string& data_file, const string& output_path, const mixMPI *mpidata) {
  chrono::time_point<chrono::high_resolution_clock> t_start = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed;
  string message;
  string timing_name;
  FILE *timing_file;
  Logger *time_logger;
  if (mpidata->rank == 0) {
    timing_name = output_path + "/c_timing_mpi"+ to_string(mpidata->rank) +".log";
    timing_file = fopen(timing_name.c_str(), "w");
    time_logger = new Logger(LOG_DEBG, timing_file);
  }
  Logger *logger = new Logger(LOG_DEBG);
  int device_count = 0;
  //===========
  // Initialise MAGMA
  //===========
#ifdef USE_MAGMA
  const magma_int_t d_array_max_size = 32; // TEMPORARY: can become configurable parameter
  magma_device_t *device_array = new magma_device_t[d_array_max_size];
  magma_int_t num_devices;
  magma_getdevices(device_array, d_array_max_size, &num_devices);
  device_count = (int)num_devices;
  delete[] device_array;
  message = "DEBUG: Proc-" + to_string(mpidata->rank) + " found " + to_string(device_count) + " GPU ";
  if (device_count > 1) message += "devices.\n";
  else message += "device.\n";
  logger->log(message, LOG_DEBG);
  logger->log("INFO: Process " + to_string(mpidata->rank) + " initializes MAGMA.\n");
  magma_int_t magma_result = magma_init();
  if (magma_result != MAGMA_SUCCESS) {
    logger->err("ERROR: Process " + to_string(mpidata->rank) + " failed to initilize MAGMA.\n");
    logger->err("PROC-" + to_string(mpidata->rank) + ": MAGMA error code " + to_string(magma_result) + "\n");
    if (mpidata->rank == 0) {
      fclose(timing_file);
      delete time_logger;
    }
    delete logger;
    return;
  }
#endif // end MAGMA initialisation
  
  //===========================
  // the following only happens on MPI process 0
  //===========================
  if (mpidata->rank == 0) {
#ifdef USE_NVTX
    nvtxRangePush("Set up");
#endif
    //=======================
    // Initialise sconf from configuration file
    //=======================
    logger->log("INFO: making legacy configuration...", LOG_INFO);
    ScattererConfiguration *sconf = NULL;
    try {
      sconf = ScattererConfiguration::from_dedfb(config_file);
    } catch(const OpenConfigurationFileException &ex) {
      logger->err("\nERROR: failed to open scatterer configuration file.\n");
      string message = "FILE: " + string(ex.what()) + "\n";
      logger->err(message);
      fclose(timing_file);
      delete time_logger;
      delete logger;
      return;
    }
    sconf->write_formatted(output_path + "/c_OEDFB");
    sconf->write_binary(output_path + "/c_TEDF");
    sconf->write_binary(output_path + "/c_TEDF.hd5", "HDF5");
    // end scatterer initialisation

    //========================
    // Initialise gconf from configuration files
    //========================
    GeometryConfiguration *gconf = NULL;
    try {
      gconf = GeometryConfiguration::from_legacy(data_file);
    } catch (const OpenConfigurationFileException &ex) {
      logger->err("\nERROR: failed to open geometry configuration file.\n");
      string message = "FILE: " + string(ex.what()) + "\n";
      logger->err(message);
      if (sconf) delete sconf;
      fclose(timing_file);
      delete time_logger;
      delete logger;
      return;
    }
    logger->log(" done.\n", LOG_INFO);
    //end gconf initialisation

#ifdef USE_NVTX
    nvtxRangePop();
#endif
    int s_nsph = sconf->number_of_spheres;
    int nsph = gconf->number_of_spheres;
    // Sanity check on number of sphere consistency, should always be verified
    if (s_nsph == nsph) {
      // Shortcuts to variables stored in configuration objects
      ScatteringAngles *p_scattering_angles = new ScatteringAngles(gconf);
      double wp = sconf->wp;
      // Open empty virtual ascii file for output
      InclusionOutputInfo *p_output = new InclusionOutputInfo(sconf, gconf, mpidata);
      InclusionIterationData *cid = new InclusionIterationData(gconf, sconf, mpidata, device_count);
      logger->log("INFO: Size of matrices to invert: " + to_string((int64_t)cid->c1->ndm) + " x " + to_string((int64_t)cid->c1->ndm) +".\n");
      time_logger->log("INFO: Size of matrices to invert: " + to_string((int64_t)cid->c1->ndm) + " x " + to_string((int64_t)cid->c1->ndm) +".\n");
      
      instr(sconf->_rcf, cid->c1);
      thdps(cid->c1->lm, cid->zpv);
      double exdc = sconf->exdc;
      double exri = sqrt(exdc);

      // Create an empty bynary file
      VirtualBinaryFile *vtppoanp = new VirtualBinaryFile();
      string tppoan_name = output_path + "/c_TPPOAN";
#ifdef USE_MAGMA
      logger->log("INFO: using MAGMA calls.\n", LOG_INFO);
#elif defined USE_LAPACK
      logger->log("INFO: using LAPACK calls.\n", LOG_INFO);
#else
      logger->log("INFO: using fall-back lucin() calls.\n", LOG_INFO);
#endif
      int iavm = gconf->iavm;
      int isam = gconf->isam;
      int inpol = gconf->in_pol;
      int nxi = sconf->number_of_scales;
      int nth = p_scattering_angles->nth;
      int nths = p_scattering_angles->nths;
      int nph = p_scattering_angles->nph;
      int nphs = p_scattering_angles->nphs;
      
      //========================
      // write a block of info to virtual binary file
      //========================
      vtppoanp->append_line(VirtualBinaryLine(iavm));
      vtppoanp->append_line(VirtualBinaryLine(isam));
      vtppoanp->append_line(VirtualBinaryLine(inpol));
      vtppoanp->append_line(VirtualBinaryLine(nxi));
      vtppoanp->append_line(VirtualBinaryLine(nth));
      vtppoanp->append_line(VirtualBinaryLine(nph));
      vtppoanp->append_line(VirtualBinaryLine(nths));
      vtppoanp->append_line(VirtualBinaryLine(nphs));
      if (sconf->idfc < 0) {
	cid->vk = cid->xip * cid->wn;
	p_output->vec_vk[0] = cid->vk;
      }
      
      // do the first iteration on jxi488 separately, since it seems to be different from the others
      int jxi488;
      int initialmaxrefiters = cid->maxrefiters;
      
//       chrono::time_point<chrono::high_resolution_clock> start_iter_1 = chrono::high_resolution_clock::now();
// #ifdef USE_NVTX
//       nvtxRangePush("First iteration");
// #endif
      // use these pragmas, which should have no effect on parallelism, just to push OMP nested levels at the same level also in the first wavelength iteration
      int jer = 0;
// #pragma omp parallel
//       {
// #pragma omp single
// 	{
// 	  jer = inclusion_jxi488_cycle(jxi488, sconf, gconf, p_scattering_angles, cid, p_output, output_path, vtppoanp);
// 	}
//       }
// #ifdef USE_NVTX
//       nvtxRangePop();
// #endif
//       chrono::time_point<chrono::high_resolution_clock> end_iter_1 = chrono::high_resolution_clock::now();
//       elapsed = start_iter_1 - t_start;
//       string message = "INFO: Calculation setup took " + to_string(elapsed.count()) + "s.\n";
//       logger->log(message);
//       time_logger->log(message);
//       elapsed = end_iter_1 - start_iter_1;
//       message = "INFO: First iteration took " + to_string(elapsed.count()) + "s.\n";
//       logger->log(message);
//       time_logger->log(message);
      /* for the next iterations, just always do maxiter iterations, assuming the accuracy is good enough */
      // cid->refinemode = 0;
      // /* add an extra iteration for margin, if this does not exceed initialmaxrefiters */
      // // if (cid->maxrefiters < initialmaxrefiters) cid->maxrefiters++;
      // if (jer != 0) {
      // 	// First loop failed. Halt the calculation.
      // 	fclose(timing_file);
      // 	delete time_logger;
      // 	delete p_output;
      // 	delete p_scattering_angles;
      // 	delete cid;
      // 	delete logger;
      // 	delete sconf;
      // 	delete gconf;
      // 	return;
      // }

      //==================================================
      // do the first outputs here, so that I open here the new files, afterwards I only append
      //==================================================
      vtppoanp->write_to_disk(output_path + "/c_TPPOAN");
      delete vtppoanp;
      
      // here go the calls that send data to be duplicated on other MPI processes from process 0 to others, using MPI broadcasts, but only if MPI is actually used
#ifdef MPI_VERSION
      if (mpidata->mpirunning) {
	gconf->mpibcast(mpidata);
	sconf->mpibcast(mpidata);	    
	cid->mpibcast(mpidata);
	p_scattering_angles->mpibcast(mpidata);
      }	
#endif
      // Create this variable and initialise it with a default here, so that it is defined anyway, with or without OpenMP support enabled
      int ompnumthreads = 1;
      // this is for MPI process 0 (or even if we are not using MPI at all)
      int myjxi488startoffset = 0;
      int myMPIstride = ompnumthreads;
      int myMPIblock = ompnumthreads;
      // Define here shared arrays of virtual ascii and binary files, so that thread 0 will be able to access them all later
      InclusionOutputInfo **p_outarray = NULL;
      VirtualBinaryFile **vtppoanarray = NULL;

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
	  p_outarray = new InclusionOutputInfo*[ompnumthreads];
	  vtppoanarray = new VirtualBinaryFile*[ompnumthreads];
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
	InclusionIterationData *cid_2 = NULL;
	InclusionOutputInfo *p_output_2 = NULL;
	VirtualBinaryFile *vtppoanp_2 = NULL;
	// for threads other than the 0, create distinct copies of all relevant data, while for thread 0 just define new references / pointers to the original ones
	if (myompthread == 0) {
	  cid_2 = cid;
	  // OMP thread 0 of MPI process 0 holds the pointer to the full output structure
	  p_output_2 = p_output;
	  p_outarray[0] = p_output_2;
	} else {
	  // this is not thread 0, so do create fresh copies of all local variables
	  cid_2 = new InclusionIterationData(*cid);
	}
	// make sure all threads align here: I don't want the following loop to accidentally start for thread 0, possibly modifying some variables before they are copied by all other threads
	if (myompthread==0) {
	  logger->log("Syncing OpenMP threads and starting the loop on wavelengths\n");
	}
#pragma omp barrier
	// ok, now I can actually start the parallel calculations
	for (int ixi488=1; ixi488<=cid_2->number_of_scales; ixi488 +=myMPIstride) {
	  // the parallel loop over MPI processes covers a different set of indices for each thread
#pragma omp barrier
	  int myjxi488 = ixi488+myompthread;
	  // each thread opens new virtual files and stores their pointers in the shared array
	  vtppoanp_2 = new VirtualBinaryFile();
	  // each thread puts a copy of the pointers to its virtual files in the shared arrays
	  vtppoanarray[myompthread] = vtppoanp_2;
#pragma omp barrier

	  // each MPI process handles a number of contiguous scales corresponding to its number of OMP threads at this omp level of parallelism
	  if (myjxi488 <= cid_2->number_of_scales) {
	    if (myompthread > 0) {
	      // UPDATE: non-0 threads need to allocate memory for one scale at a time.
	      p_output_2 = new InclusionOutputInfo(sconf, gconf, mpidata, myjxi488, 1);
	      p_outarray[myompthread] = p_output_2;
	    }
	    int jer = inclusion_jxi488_cycle(myjxi488, sconf, gconf, p_scattering_angles, cid_2, p_output_2, output_path, vtppoanp_2);
	  } else {
	    if (myompthread > 0) {
	      // If there is no input for this thread, set output pointer to NULL.
	      p_outarray[myompthread] = new InclusionOutputInfo(1);
	    }
	  }
#pragma omp barrier

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
	      vtppoanarray[0]->append(*(vtppoanarray[ti]));
	      delete vtppoanarray[ti];
	    }
	  }
#pragma omp barrier
	  //==============================================
	  // Collect all virtual files on thread 0 of MPI process 0, and append them to disk
	  //==============================================
	  if (myompthread == 0) {
	    // thread 0 writes its virtual files, now including contributions from all threads, to disk, and deletes them
	    // p_outarray[0]->append_to_disk(output_path + "/c_OINCLU");
	    // delete p_outarray[0];
	    vtppoanarray[0]->append_to_disk(output_path + "/c_TPPOAN");
	    delete vtppoanarray[0];

#ifdef MPI_VERSION
	    if (mpidata->mpirunning) {
	      // only go through this if MPI has been actually used
	      for (int rr=1; rr<mpidata->nprocs; rr++) {
		// get the data from process rr by receiving it in total memory structure
		p_outarray[0]->mpireceive(mpidata, rr);
		// get the data from process rr, creating a new virtual ascii file
		// VirtualAsciiFile *p_output = new VirtualAsciiFile(mpidata, rr);
		// append to disk and delete virtual ascii file
		// p_output->append_to_disk(output_path + "/c_OINCLU");
		// delete p_output;
		
		// get the data from process rr, creating a new virtual binary file
		VirtualBinaryFile *vtppoanp = new VirtualBinaryFile(mpidata, rr);
		// append to disk and delete virtual binary file
		vtppoanp->append_to_disk(output_path + "/c_TPPOAN");
		delete vtppoanp;
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
	  delete[] p_outarray;
	  delete[] vtppoanarray;
	}
	{
	  string message = "INFO: Closing thread-local output files of thread " + to_string(myompthread) + " and syncing threads.\n";
	  logger->log(message);
	}
#ifdef USE_NVTX
	nvtxRangePop();
#endif
	delete cid_2;
      }
      delete p_scattering_angles;
      p_output->write(output_path + "/c_OINCLU.hd5", "HDF5");
      p_output->write(output_path + "/c_OINCLU", "LEGACY");
      delete p_output;
    } // closes s_nsph == nsph check
	  
    else { // Sphere number inconsistency error.
      throw UnrecognizedConfigurationException(
	"Inconsistent geometry and scatterer configurations."
      );
    }

    delete sconf;
    delete gconf;
#ifdef USE_MAGMA
    logger->log("INFO: Process " + to_string(mpidata->rank) + " finalizes MAGMA.\n");
    magma_finalize();
#endif
    chrono::time_point<chrono::high_resolution_clock> t_end = chrono::high_resolution_clock::now();
    elapsed = t_end - t_start;
    string message = "INFO: Calculation lasted " + to_string(elapsed.count()) + "s.\n";
    logger->log(message);
    logger->log("Finished: output written to " + output_path + "/c_OINCLU\n");
    time_logger->log(message);
    fclose(timing_file);
    delete time_logger;
    delete logger;
  } // end instructions block of MPI process 0
  
    //===============================
    // instruction block for MPI processes different from 0
    //===============================
#ifdef MPI_VERSION
  else {
    // here go the code for MPI processes other than 0
    // copy gconf, sconf, cid and p_scattering_angles from MPI process 0
    GeometryConfiguration *gconf = new GeometryConfiguration(mpidata);
    ScattererConfiguration *sconf = new ScattererConfiguration(mpidata);
    InclusionIterationData *cid = new InclusionIterationData(mpidata, device_count);
    ScatteringAngles *p_scattering_angles = new ScatteringAngles(mpidata);

    // Create this variable and initialise it with a default here, so that it is defined anyway, with or without OpenMP support enabled
    int ompnumthreads = 1;
    InclusionOutputInfo **p_outarray = NULL;
    VirtualBinaryFile **vtppoanarray = NULL;
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
	p_outarray = new InclusionOutputInfo*[ompnumthreads];
	vtppoanarray = new VirtualBinaryFile*[ompnumthreads];
      }
#pragma omp barrier
      // To test parallelism, I will now start feeding this function with "clean" copies of the parameters, so that they will not be changed by previous iterations, and each one will behave as the first one. Define all (empty) variables here, so they have the correct scope, then they get different definitions depending on thread number
      InclusionIterationData *cid_2 = NULL;
      InclusionOutputInfo *p_output_2 = NULL;
      VirtualBinaryFile *vtppoanp_2 = NULL;
      // PLACEHOLDER
      // for threads other than the 0, create distinct copies of all relevant data, while for thread 0 just define new references / pointers to the original ones
      if (myompthread == 0) {
	cid_2 = cid;
      } else {
	// this is not thread 0, so do create fresh copies of all local variables
	cid_2 = new InclusionIterationData(*cid);
      }
      // make sure all threads align here: I don't want the following loop to accidentally start for thread 0, possibly modifying some variables before they are copied by all other threads
#pragma omp barrier
      // ok, now I can actually start the parallel calculations
      for (int ixi488=1; ixi488<=cid_2->number_of_scales; ixi488 +=myMPIstride) {
	// the parallel loop over MPI processes covers a different set of indices for each thread
#pragma omp barrier
	int myjxi488 = ixi488 + myjxi488startoffset + myompthread;
	// each thread opens new virtual files and stores their pointers in the shared array
	vtppoanp_2 = new VirtualBinaryFile();
	// each thread puts a copy of the pointers to its virtual files in the shared arrays
	vtppoanarray[myompthread] = vtppoanp_2;
#pragma omp barrier
	if (myompthread==0) logger->log("Syncing OpenMP threads and starting the loop on wavelengths\n");
	// ok, now I can actually start the parallel calculations
	// each MPI process handles a number of contiguous scales corresponding to its number of OMP threads at this omp level of parallelism
	if (myjxi488 <= cid_2->number_of_scales) {
	  if (myompthread > 0) {
	    // UPDATE: non-0 threads need to allocate memory for one scale at a time.
	    p_output_2 = new InclusionOutputInfo(sconf, gconf, mpidata, myjxi488, 1);
	    p_outarray[myompthread] = p_output_2;
	  } else {
	    // Thread 0 of non-zero MPI processes needs to allocate memory for the
	    // output of all threads _doing something_.
	    int iterstodo = cid_2->number_of_scales - myjxi488 + 1;
	    if (iterstodo > ompnumthreads) iterstodo = ompnumthreads;
	    p_output_2 = new InclusionOutputInfo(sconf, gconf, mpidata, myjxi488, iterstodo);
	    p_outarray[0] = p_output_2;
	  }
	  int jer = inclusion_jxi488_cycle(myjxi488, sconf, gconf, p_scattering_angles, cid_2, p_output_2, output_path, vtppoanp_2);
	} else {
	  p_outarray[myompthread] = new InclusionOutputInfo(1);
	}

#pragma omp barrier
	// threads different from 0 append their virtual files to the one of thread 0, and delete them
	if (myompthread == 0) {
	  for (int ti=1; ti<ompnumthreads; ti++) {
	    p_outarray[0]->insert(*(p_outarray[ti]));
	    delete p_outarray[ti];
	    p_outarray[ti] = NULL;
	    vtppoanarray[0]->append(*(vtppoanarray[ti]));
	    delete vtppoanarray[ti];
	  }
	  // thread 0 sends the collected virtualfiles to thread 0 of MPI process 0, then deletes them
	  for (int rr=1; rr<mpidata->nprocs; rr++) {
	    if (rr == mpidata->rank) {
	      p_outarray[0]->mpisend(mpidata);
	      delete p_outarray[0];
	      p_outarray[0] = NULL;
	      vtppoanarray[0]->mpisend(mpidata);
	      delete vtppoanarray[0];
	    }
	    int test = MPI_Barrier(MPI_COMM_WORLD);
	  }
	}
      } // close strided loop running on MPI processes
      
	// Clean memory
#pragma omp barrier
      if (myompthread == 0) {
	delete[] p_outarray;
	delete[] vtppoanarray;
      }
      delete cid_2;

    } // close pragma omp parallel
    delete p_scattering_angles;
    delete sconf;
    delete gconf;
#endif
#ifdef USE_MAGMA
    logger->log("INFO: Process " + to_string(mpidata->rank) + " finalizes MAGMA.\n");
    magma_finalize();
#endif
    delete logger;
#ifdef MPI_VERSION
  }
#endif
}

int inclusion_jxi488_cycle(int jxi488, ScattererConfiguration *sconf, GeometryConfiguration *gconf, ScatteringAngles *sa, InclusionIterationData *cid, InclusionOutputInfo *output, const string& output_path, VirtualBinaryFile *vtppoanp) {
  int nxi = sconf->number_of_scales;
  const dcomplex cc0 = 0.0 + I * 0.0;
  const double pi = acos(-1.0);
  char virtual_line[256];
  string message = "INFO: running scale iteration " + to_string(jxi488) + " of " + to_string(nxi) + ".\n";
  Logger *logger = new Logger(LOG_DEBG);
  logger->log(message);
  chrono::duration<double> elapsed;
  chrono::time_point<chrono::high_resolution_clock> interval_start, interval_end;
  int jer = 0;
  int lcalc = 0;
  int jaw = 1;
  // this is now part of cid, so don't mess with it here, just copy it by reference
  bool &is_first_scale = cid->is_first_scale;
  //bool is_first_scale = (jxi488 == 1);
  // int li = cid->c1->li;
  // int le = cid->c1->le;
  // int lm = cid->c1->lm;
  const int nsph = cid->c1->nsph;
  np_int mxndm = gconf->mxndm;
  const int iavm = gconf->iavm;
  const int inpol = gconf->in_pol;
  const int npnt = cid->c1->npnt;
  const int npntts = cid->c1->npntts;
  const int isam = gconf->iavm;
  const int jwtm = gconf->jwtm;
  int isq, ibf;
  int last_configuration;
  dcomplex ent, entn;
  double enti;
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
  int idfc = (int)sconf->idfc;
  double vkarg = 0.0;
  if (idfc >= 0) {
    cid->vk = xi * cid->wn;
    vkarg = cid->vk;
    output->vec_vk[jindex - 1] = cid->vk;
    output->vec_xi[jindex - 1] = xi;
    // goes to 120
  } else { // label 119
    vkarg = xi * cid->vk;
    cid->sqsfi = 1.0 / (xi * xi);
    output->vec_vk[jindex - 1] = cid->vk;
    output->vec_xi[jindex - 1] = xi;
  }
  // Dynamic order check
  const int max_li = gconf->li;
  const int max_le = gconf->le;
  const double alamb = 2.0 * pi / cid->vk;
  double size_par_li = 2.0 * pi * sqrt(exdc) * sconf->get_max_radius() / alamb;
  int recommended_li = 4 + (int)ceil(size_par_li + 4.05 * pow(size_par_li, 1.0 / 3.0));
  double size_par_le = 2.0 * pi * sqrt(exdc) * sconf->get_particle_radius(gconf) / alamb;
  int recommended_le = 1 + (int)ceil(size_par_le + 11.0 * pow(size_par_le, 1.0 / 3.0));
  if (recommended_li != cid->c1->li || recommended_le != cid->c1->le) {
    if (recommended_li > cid->c1->li) {
      message = "WARNING: internal order " + to_string(cid->c1->li) + " for scale iteration "
	+ to_string(jxi488) + " too low (recommended order is " + to_string(recommended_li)
	+ ").\n";
      logger->log(message, LOG_WARN);
    } else if (recommended_li < cid->c1->li) {
      if (gconf->dyn_order_flag > 0) {
	message = "INFO: lowering internal order from " + to_string(cid->c1->li) + " to "
	  + to_string(recommended_li) + " for scale iteration " + to_string(jxi488) + ".\n";
	logger->log(message, LOG_INFO);
      } else {
	message = "WARNING: internal order " + to_string(cid->c1->li) + " too high for scale iteration "
	  + to_string(jxi488) + ".\n";
	logger->log(message, LOG_WARN);
      }
    }
    if (recommended_le > cid->c1->le) {
      message = "WARNING: external order " + to_string(cid->c1->le) + " for scale iteration "
	+ to_string(jxi488) + " too low (recommended order is " + to_string(recommended_le)
	+ ").\n";
      logger->log(message, LOG_WARN);
    } else if (recommended_le < cid->c1->le) {
      if (gconf->dyn_order_flag > 0) {
	message = "INFO: lowering external order from " + to_string(cid->c1->le) + " to "
	  + to_string(recommended_le) + " for scale iteration " + to_string(jxi488) + ".\n";
	logger->log(message, LOG_INFO);
      } else {
	message = "WARNING: external order " + to_string(cid->c1->le) + " too high for scale iteration "
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
	cid->refinemode = 2;
      }
    }
  }
  int li = cid->c1->li;
  int le = cid->c1->le;
  int lm = cid->c1->lm;
  // End of dynamic order check
  // label 120
  double sze = vkarg * cid->extr;
  last_configuration = 0;
  for (int i133 = 1; i133 <= cid->c1->nsph; i133++) {
    int iogi = cid->c1->iog[i133 - 1];
    if (iogi != i133) {
      for (int l123 = 1; l123 <= cid->c1->li; l123++) {
	cid->c1->rmi[l123 - 1][i133 - 1] = cid->c1->rmi[l123 - 1][iogi - 1];
	cid->c1->rei[l123 - 1][i133 - 1] = cid->c1->rei[l123 - 1][iogi - 1];
      } // l123 loop
    } else { // label 125
      last_configuration++;
      int nsh = cid->c1->nshl[last_configuration - 1];
      int ici = (nsh + 1) / 2;
      if (i133 == 1) ici++;
      if (idfc == 0) {
	for (int ic = 0; ic < ici; ic++)
	  cid->c1->dc0[ic] = sconf->get_dielectric_constant(ic, i133 - 1, jxi488 - 1);
	// goes to 129
      } else { // label 127
	if (is_first_scale) {
	  for (int ic = 0; ic < ici; ic++) {
	    cid->c1->dc0[ic] = sconf->get_dielectric_constant(ic, i133 - 1, 0);
	  }
	}
      }
      // label 129
      if (i133 == 1) {
	ent = cid->c1->dc0[ici - 1];
	enti = imag(ent);
	entn = csqrt(ent);
	// goes to 131
      } else { // label 130
	if (nsh % 2 == 0) cid->c1->dc0[ici] = ent;
      }
      indme(i133, npnt, npntts, vkarg, ent, enti, entn, jer, lcalc, cid->arg, cid->c1);
      if (jer != 0) {
	output->vec_ier[jindex - 1] = 1;
	output->vec_jxi[jindex - 1] = -jxi488;
	message = "ERROR: indme failed with error code " + to_string(jer) + ".\n";
	logger->log(message, LOG_ERRO);
	delete logger;
	return jer;
	//break;
      }
    }
  } // i133 loop
  ospv(cid->c1, vkarg, sze, exri, entn, enti, jer, lcalc, cid->arg);
  if (jer != 0) {
    output->vec_ier[jindex - 1] = 2;
    output->vec_jxi[jindex - 1] = -jxi488;
    message = "ERROR: ospv failed with error code " + to_string(jer) + ".\n";
    logger->log(message, LOG_ERRO);
    delete logger;
    return jer;
    // break;
  } // i133 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  interval_start = chrono::high_resolution_clock::now();
#ifdef USE_NVTX
  nvtxRangePush("Calculate inverted matrix");
#endif
#ifdef DEBUG_AM
  /* now, before cms, output am to p_outam0 */
  VirtualAsciiFile *outam0 = new VirtualAsciiFile();
  string outam0_name = output_path + "/c_AM0_JXI" + to_string(jxi488) + ".txt";
  sprintf(virtual_line, " AM matrix before CMS\n");
  outam0->append_line(virtual_line);
  sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
  outam0->append_line(virtual_line);
  write_dcomplex_matrix(outam0, cid->am, cid->c1->ndm, cid->c1->ndm);
  outam0->write_to_disk(outam0_name);
  delete outam0;
#endif
  incms(cid->am, enti, cid->c1);
#ifdef DEBUG_AM
  VirtualAsciiFile *outam1 = new VirtualAsciiFile();
  string outam1_name = output_path + "/c_AM1_JXI" + to_string(jxi488) + ".txt";
  sprintf(virtual_line, " AM matrix after CMS before LUCIN\n");
  outam1->append_line(virtual_line);
  sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
  outam1->append_line(virtual_line);
  write_dcomplex_matrix(outam1, cid->am, cid->c1->ndm, cid->c1->ndm, " %5d %5d (%17.8lE,%17.8lE)\n", 1);
  outam1->write_to_disk(outam1_name);
  delete outam1;
#endif
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  interval_end = chrono::high_resolution_clock::now();
  elapsed = interval_end - interval_start;
  message = "INFO: matrix calculation for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
  logger->log(message);
  interval_start = chrono::high_resolution_clock::now();
#ifdef USE_NVTX
  nvtxRangePush("Invert the matrix");
#endif
  // we the accuracygoal in, get the actual accuracy back out
  double actualaccuracy = cid->accuracygoal;
  invert_matrix(cid->am, cid->c1->ndm, jer, cid->maxrefiters, actualaccuracy, cid->refinemode, output_path, jxi488, mxndm, cid->proc_device);
  // in principle, we should check whether the returned actualaccuracy is indeed lower than the accuracygoal, and do something about it if not
  if (gconf->refine_flag > 0) {
    if (cid->refinemode==2) {
      message = "DEBUG: iterative refinement enabled at run-time.\n";
      logger->log(message, LOG_DEBG);
      message = "INFO: calibration obtained accuracy " + to_string(actualaccuracy) + " (" + to_string(cid->accuracygoal) + " requested) in " + to_string(cid->maxrefiters) + " refinement iterations\n";
      logger->log(message);
      if (actualaccuracy > 1e-2) {
	printf("Accuracy worse than 0.01, stopping");
	exit(1);
      }
    }
  }
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  interval_end = chrono::high_resolution_clock::now();
  elapsed = interval_end - interval_start;
  message = "INFO: matrix inversion for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
  logger->log(message);
  if (jer != 0) {
    message = "ERROR: matrix inversion ended with error code " + to_string(jer) + ".\n";
    logger->err(message);
    delete logger;
    return jer;
    // break; // jxi488 loop: goes to memory clean
  }
  interval_start = chrono::high_resolution_clock::now();
#ifdef USE_NVTX
  nvtxRangePush("Average calculation");
#endif
  exma(cid->am, cid->c1);
#ifdef DEBUG_AM
  VirtualAsciiFile *outam3 = new VirtualAsciiFile();
  string outam3_name = output_path + "/c_AM3_JXI" + to_string(jxi488) + ".txt";
  sprintf(virtual_line, " AM matrix after EXMA\n");
  outam3->append_line(virtual_line);
  sprintf(virtual_line, " I1+1   I2+1    Real    Imag\n");
  outam3->append_line(virtual_line);
  write_dcomplex_matrix(outam3, cid->am, cid->c1->ndm, cid->c1->ndm);
  outam3->write_to_disk(outam3_name);
  delete outam3;
#endif
  if (idfc >= 0) {
    if (jxi488 == jwtm) {
      int nlemt = 2 * cid->c1->nlem;
      string ttms_name = output_path + "/c_TTMS.hd5";
      TransitionMatrix::write_binary(ttms_name, nlemt, lm, cid->vk, exri, cid->c1->am0m, "HDF5");
      ttms_name = output_path + "/c_TTMS";
      TransitionMatrix::write_binary(ttms_name, nlemt, lm, cid->vk, exri, cid->c1->am0m);
    }
  }
  // label 156: continue from here
  last_configuration = 0;
  for (int i168 = 1; i168 <= nsph; i168++) {
    if (cid->c1->iog[i168 - 1] >= i168) {
      int i = i168 - 1;
      last_configuration++;
      oindex = (jindex - 1) * (num_configs + 1) + last_configuration - 1;
      if (cid->c1->nshl[i168 - 1] != 1) {
	output->vec_sphere_ref_indices[oindex] = cc0;
	output->vec_sphere_sizes[oindex] = cid->c1->vsz[i];
      } else {
	output->vec_sphere_ref_indices[oindex] = cid->c1->vkt[i];
	output->vec_sphere_sizes[oindex] = cid->c1->vsz[i];
      }
    } 
  } // i168 loop
  oindex = (jindex - 1) * (num_configs + 1) + num_configs;
  output->vec_sphere_sizes[oindex] = sze;
  output->vec_sphere_ref_indices[oindex] = entn;
  // label 160
  double cs0 = 0.25 * cid->vk * cid->vk * cid->vk / acos(0.0);
  double csch = 2.0 * cid->vk * cid->sqsfi / cid->c1->gcs;
  double sqk = cid->vk * cid->vk * exdc;
  vtppoanp->append_line(VirtualBinaryLine(cid->vk));
  pcrsm0(cid->vk, exri, inpol, cid->c1);
  apcra(cid->zpv, cid->c1->le, cid->c1->am0m, inpol, sqk, cid->gapm, cid->gappm);
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  interval_end = chrono::high_resolution_clock::now();
  elapsed = interval_end - interval_start;
  message = "INFO: average calculation for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
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
      } else { // label 180, NK == 1 AND JXI488 == 1
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
	if (isam > 1) ths += sa->thsca;
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
	    phs = sa->ph + phsph;
	    if (phs > 360.0) phs -= 360.0;
	  }
	  // label 188
	  bool goto190 = (sa->nks == 1 && (!(is_first_scale) || jth486 > 1 || jph484 > 1));
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
	  vtppoanp->append_line(VirtualBinaryLine(th));
	  vtppoanp->append_line(VirtualBinaryLine(ph));
	  vtppoanp->append_line(VirtualBinaryLine(ths));
	  vtppoanp->append_line(VirtualBinaryLine(phs));
	  vtppoanp->append_line(VirtualBinaryLine(cid->scan));
	  if (jaw != 0) {
	    jaw = 0;
	    mextc(cid->vk, exri, cid->c1->fsacm, cid->cextlr, cid->cext);
	    // We now have some implicit loops writing to binary
	    for (int i = 0; i < 4; i++) {
	      for (int j = 0; j < 4; j++) {
		double value = cid->cext[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    for (int i = 0; i < 2; i++) {
	      double value = cid->c1->scscm[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = real(cid->c1->scscpm[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = imag(cid->c1->scscpm[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = cid->c1->ecscm[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = real(cid->c1->ecscpm[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = imag(cid->c1->ecscpm[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	    }
	    for (int i = 0; i < 3; i++) {
	      for (int j = 0; j < 2; j++) {
		double value = cid->gapm[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = real(cid->gappm[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = imag(cid->gappm[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    int jlr = 2;
	    for (int ilr210 = 1; ilr210 <= 2; ilr210++) {
	      int ipol = (ilr210 % 2 == 0) ? 1 : -1;
	      if (ilr210 == 2) jlr = 1;
	      double extsm = cid->c1->ecscm[ilr210 - 1];
	      double qextm = extsm * cid->sqsfi / cid->c1->gcs;
	      double scasm = cid->c1->scscm[ilr210 - 1];
	      double albdm = scasm / extsm;
	      double qscam = scasm * cid->sqsfi / cid->c1->gcs;
	      double abssm = extsm - scasm;
	      double qabsm = abssm * cid->sqsfi / cid->c1->gcs;
	      dcomplex s0m = cid->c1->fsacm[ilr210 - 1][ilr210 - 1] * exri;
	      double qschum = imag(s0m) * csch;
	      double pschum = real(s0m) * csch;
	      double s0magm = cabs(s0m) * cs0;
	      // if (inpol == 0) {
	      // sprintf(virtual_line, "   LIN %2d\n", ipol);
	      // output->append_line(virtual_line);
	      // } else { // label 206
	      // sprintf(virtual_line, "  CIRC %2d\n", ipol);
	      // output->append_line(virtual_line);
	      // }
	      // label 208
	      double rapr = cid->c1->ecscm[ilr210 - 1] - cid->gapm[2][ilr210 - 1];
	      double cosav = cid->gapm[2][ilr210 - 1] / cid->c1->scscm[ilr210 - 1];
	      double fz = rapr;
	      if (ipol == -1) {
		output->vec_scs1[jindex - 1] = scasm;
		output->vec_abs1[jindex - 1] = abssm;
		output->vec_exs1[jindex - 1] = extsm;
		output->vec_albeds1[jindex - 1] = albdm;
		output->vec_scsrt1[jindex - 1] = qscam;
		output->vec_absrt1[jindex - 1] = qabsm;
		output->vec_exsrt1[jindex - 1] = qextm;
		output->vec_fsas11[jindex - 1] = cid->c1->fsacm[0][0];
		output->vec_fsas21[jindex - 1] = cid->c1->fsacm[1][0];
		output->vec_qschu1[jindex -1] = qschum;
		output->vec_pschu1[jindex -1] = pschum;
		output->vec_s0mag1[jindex -1] = s0magm;
		output->vec_cosav1[jindex -1] = cosav;
		output->vec_raprs1[jindex -1] = rapr;
		output->vec_fk1[jindex -1] = fz;
	      } else if (ipol == 1) {
		output->vec_scs2[jindex - 1] = scasm;
		output->vec_abs2[jindex - 1] = abssm;
		output->vec_exs2[jindex - 1] = extsm;
		output->vec_albeds2[jindex - 1] = albdm;
		output->vec_scsrt2[jindex - 1] = qscam;
		output->vec_absrt2[jindex - 1] = qabsm;
		output->vec_exsrt2[jindex - 1] = qextm;
		output->vec_fsas22[jindex - 1] = cid->c1->fsacm[1][1];
		output->vec_fsas12[jindex - 1] = cid->c1->fsacm[0][1];
		output->vec_qschu2[jindex -1] = qschum;
		output->vec_pschu2[jindex -1] = pschum;
		output->vec_s0mag2[jindex -1] = s0magm;
		output->vec_cosav2[jindex -1] = cosav;
		output->vec_raprs2[jindex -1] = rapr;
		output->vec_fk2[jindex -1] = fz;
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
	  // label 220
	  pcros(cid->vk, exri, cid->c1);
	  mextc(cid->vk, exri, cid->c1->fsac, cid->cextlr, cid->cext);
	  mmulc(cid->c1->vint, cid->cmullr, cid->cmul);
	  if (jw != 0) {
	    jw = 0;
	    // Some implicit loops writing to binary.
	    for (int i = 0; i < 4; i++) {
	      for (int j = 0; j < 4; j++) {
		double value = cid->cext[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    for (int i = 0; i < 2; i++) {
	      double value = cid->c1->scsc[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = real(cid->c1->scscp[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = imag(cid->c1->scscp[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = cid->c1->ecsc[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = real(cid->c1->ecscp[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = imag(cid->c1->ecscp[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	    }
	    for (int i = 0; i < 3; i++) {
	      for (int j = 0; j < 2; j++) {
		double value = cid->gap[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = real(cid->gapp[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = imag(cid->gapp[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    for (int i = 0; i < 2; i++) {
	      for (int j = 0; j < 3; j++) {
		double value = cid->tqce[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = real(cid->tqcpe[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = imag(cid->tqcpe[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    for (int i = 0; i < 2; i++) {
	      for (int j = 0; j < 3; j++) {
		double value = cid->tqcs[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = real(cid->tqcps[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
		value = imag(cid->tqcps[i][j]);
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    for (int i = 0; i < 3; i++) {
	      double value = cid->u[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = cid->up[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = cid->un[i];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	    }
	  }
	  // label 254
	  for (int i = 0; i < 16; i++) {
	    double value = real(cid->c1->vint[i]);
	    vtppoanp->append_line(VirtualBinaryLine(value));
	    value = imag(cid->c1->vint[i]);
	    vtppoanp->append_line(VirtualBinaryLine(value));
	  }
	  for (int i = 0; i < 4; i++) {
	    for (int j = 0; j < 4; j++) {
	      double value = cid->cmul[i][j];
	      vtppoanp->append_line(VirtualBinaryLine(value));
	    }
	  }
	  oindex = (jindex - 1) * ndirs + done_dirs;
	  int jlr = 2;
	  for (int ilr290 = 1; ilr290 <= 2; ilr290++) {
	    int ipol = (ilr290 % 2 == 0) ? 1 : -1;
	    if (ilr290 == 2) jlr = 1;
	    double extsec = cid->c1->ecsc[ilr290 - 1];
	    double qext = extsec * cid->sqsfi / cid->c1->gcs;
	    double scasec = cid->c1->scsc[ilr290 - 1];
	    double albedc = scasec / extsec;
	    double qsca = scasec * cid->sqsfi / cid->c1->gcs;
	    double abssec = extsec - scasec;
	    double qabs = abssec * cid->sqsfi / cid->c1->gcs;
	    dcomplex s0 = cid->c1->fsac[ilr290 - 1][ilr290 - 1] * exri;
	    double qschu = imag(s0) * csch;
	    double pschu = real(s0) * csch;
	    double s0mag = cabs(s0) * cs0;
	    // label 275
	    if (ipol == -1) {
	      output->vec_dir_scs1[oindex] = scasec;
	      output->vec_dir_abs1[oindex] = abssec;
	      output->vec_dir_exs1[oindex] = extsec;
	      output->vec_dir_albeds1[oindex] = albedc;
	      output->vec_dir_scsrt1[oindex] = qsca;
	      output->vec_dir_absrt1[oindex] = qabs;
	      output->vec_dir_exsrt1[oindex] = qext;
	      output->vec_dir_fsas11[oindex] = cid->c1->fsac[0][0];
	      output->vec_dir_fsas21[oindex] = cid->c1->fsac[1][0];
	      output->vec_dir_sas11[oindex] = cid->c1->sac[0][0];
	      output->vec_dir_sas21[oindex] = cid->c1->sac[1][0];
	      output->vec_dir_qschu1[oindex] = qschu;
	      output->vec_dir_pschu1[oindex] = pschu;
	      output->vec_dir_s0mag1[oindex] = s0mag;
	    } else if (ipol == 1) {
	      output->vec_dir_scs2[oindex] = scasec;
	      output->vec_dir_abs2[oindex] = abssec;
	      output->vec_dir_exs2[oindex] = extsec;
	      output->vec_dir_albeds2[oindex] = albedc;
	      output->vec_dir_scsrt2[oindex] = qsca;
	      output->vec_dir_absrt2[oindex] = qabs;
	      output->vec_dir_exsrt2[oindex] = qext;
	      output->vec_dir_fsas22[oindex] = cid->c1->fsac[1][1];
	      output->vec_dir_fsas12[oindex] = cid->c1->fsac[0][1];
	      output->vec_dir_sas22[oindex] = cid->c1->sac[1][1];
	      output->vec_dir_sas12[oindex] = cid->c1->sac[0][1];
	      output->vec_dir_qschu2[oindex] = qschu;
	      output->vec_dir_pschu2[oindex] = pschu;
	      output->vec_dir_s0mag2[oindex] = s0mag;
	    }
	    bool goto290 = isam >= 0 && (jths > 1 || jphs > 1);
	    if (!goto290) {
	      cid->gapv[0] = cid->gap[0][ilr290 - 1];
	      cid->gapv[1] = cid->gap[1][ilr290 - 1];
	      cid->gapv[2] = cid->gap[2][ilr290 - 1];
	      double extins = cid->c1->ecsc[ilr290 - 1];
	      double scatts = cid->c1->scsc[ilr290 - 1];
	      double rapr, cosav, fp, fn, fk, fx, fy, fz;
	      rftr(cid->u, cid->up, cid->un, cid->gapv, extins, scatts, rapr, cosav, fp, fn, fk, fx, fy, fz);
	      if (ipol == -1) {
		output->vec_dir_cosav1[oindex] = cosav;
		output->vec_dir_rapr1[oindex] = rapr;
		output->vec_dir_fl1[oindex] = fp;
		output->vec_dir_fr1[oindex] = fn;
		output->vec_dir_fk1[oindex] = fk;
		output->vec_dir_fx1[oindex] = fx;
		output->vec_dir_fy1[oindex] = fy;
		output->vec_dir_fz1[oindex] = fz;
	      } else if (ipol == 1) {
		output->vec_dir_cosav2[oindex] = cosav;
		output->vec_dir_rapr2[oindex] = rapr;
		output->vec_dir_fl2[oindex] = fp;
		output->vec_dir_fr2[oindex] = fn;
		output->vec_dir_fk2[oindex] = fk;
		output->vec_dir_fx2[oindex] = fx;
		output->vec_dir_fy2[oindex] = fy;
		output->vec_dir_fz2[oindex] = fz;
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
		output->vec_dir_tqel1[oindex] = tep;
		output->vec_dir_tqer1[oindex] = ten;
		output->vec_dir_tqek1[oindex] = tek;
		output->vec_dir_tqsl1[oindex] = tsp;
		output->vec_dir_tqsr1[oindex] = tsn;
		output->vec_dir_tqsk1[oindex] = tsk;
		output->vec_dir_tqex1[oindex] = cid->tqce[0][0];
		output->vec_dir_tqey1[oindex] = cid->tqce[0][1];
		output->vec_dir_tqez1[oindex] = cid->tqce[0][2];
		output->vec_dir_tqsx1[oindex] = cid->tqcs[0][0];
		output->vec_dir_tqsy1[oindex] = cid->tqcs[0][1];
		output->vec_dir_tqsz1[oindex] = cid->tqcs[0][2];
	      } else if (ipol == 1) {
		output->vec_dir_tqel2[oindex] = tep;
		output->vec_dir_tqer2[oindex] = ten;
		output->vec_dir_tqek2[oindex] = tek;
		output->vec_dir_tqsl2[oindex] = tsp;
		output->vec_dir_tqsr2[oindex] = tsn;
		output->vec_dir_tqsk2[oindex] = tsk;
		output->vec_dir_tqex2[oindex] = cid->tqce[1][0];
		output->vec_dir_tqey2[oindex] = cid->tqce[1][1];
		output->vec_dir_tqez2[oindex] = cid->tqce[1][2];
		output->vec_dir_tqsx2[oindex] = cid->tqcs[1][0];
		output->vec_dir_tqsy2[oindex] = cid->tqcs[1][1];
		output->vec_dir_tqsz2[oindex] = cid->tqcs[1][2];
	      }
	    }
	  } //ilr290 loop
	  for (int i = 0; i < 4; i++) {
	    oindex = 16 * (jindex - 1) * ndirs + 16 * done_dirs + 4 * i;
	    output->vec_dir_mull[oindex] = cid->cmul[i][0];
	    output->vec_dir_mull[oindex + 1] = cid->cmul[i][1];
	    output->vec_dir_mull[oindex + 2] = cid->cmul[i][2];
	    output->vec_dir_mull[oindex + 3] = cid->cmul[i][3];
	  }
	  for (int i = 0; i < 4; i++) {
	    oindex = 16 * (jindex - 1) * ndirs + 16 * done_dirs + 4 * i;
	    output->vec_dir_mulllr[oindex] = cid->cmullr[i][0];
	    output->vec_dir_mulllr[oindex + 1] = cid->cmullr[i][1];
	    output->vec_dir_mulllr[oindex + 2] = cid->cmullr[i][2];
	    output->vec_dir_mulllr[oindex + 3] = cid->cmullr[i][3];
	  }
	  if (iavm != 0) {
	    mmulc(cid->c1->vintm, cid->cmullr, cid->cmul);
	    // Some implicit loops writing to binary.
	    for (int i = 0; i < 16; i++) {
	      double value;
	      value = real(cid->c1->vintm[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	      value = imag(cid->c1->vintm[i]);
	      vtppoanp->append_line(VirtualBinaryLine(value));
	    }
	    for (int i = 0; i < 4; i++) {
	      for (int j = 0; j < 4; j++) {
		double value = cid->cmul[i][j];
		vtppoanp->append_line(VirtualBinaryLine(value));
	      }
	    }
	    // label 318
	    for (int i = 0; i < 4; i++) {
	      oindex = 16 * (jindex - 1) + 4 * i; // if IAVM fails, try adding directions
	      output->vec_dir_mull[oindex] = cid->cmul[i][0];
	      output->vec_dir_mull[oindex + 1] = cid->cmul[i][1];
	      output->vec_dir_mull[oindex + 2] = cid->cmul[i][2];
	      output->vec_dir_mull[oindex + 3] = cid->cmul[i][3];
	    }
	    for (int i = 0; i < 4; i++) {
	      oindex = 16 * (jindex - 1) + 4 * i; // if IAVM fails, try adding directions
	      output->vec_dir_mulllr[oindex] = cid->cmul[i][0];
	      output->vec_dir_mulllr[oindex + 1] = cid->cmullr[i][1];
	      output->vec_dir_mulllr[oindex + 2] = cid->cmullr[i][2];
	      output->vec_dir_mulllr[oindex + 3] = cid->cmullr[i][3];
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
  message = "INFO: angle loop for scale " + to_string(jxi488) + " took " + to_string(elapsed.count()) + "s.\n";
  logger->log(message);
  
  logger->log("INFO: finished scale iteration " + to_string(jxi488) + " of " + to_string(nxi) + ".\n");

  delete logger;

  return jer;
}

// >>> IMPLEMENTATION OF InclusionIterationData CLASS <<< //
InclusionIterationData::InclusionIterationData(GeometryConfiguration *gconf, ScattererConfiguration *sconf, const mixMPI *mpidata, const int device_count) {
  c1 = new ParticleDescriptorInclusion(gconf, sconf);
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
  vec_zpv = new double[c1->lm * 12]();
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double**[12];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[4];
      zpv[zi][zj][0] = vec_zpv + (zi * 12) + (zj * 4);
      zpv[zi][zj][1] = vec_zpv + (zi * 12) + (zj * 4) + 2;
    }
  }
  am_vector = new dcomplex[c1->ndm * c1->ndm]();
  am = new dcomplex*[c1->ndm];
  for (int ai = 0; ai < c1->ndm; ai++) {
    am[ai] = (am_vector + ai * c1->ndm);
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

  nimd = c1->nshl[0] + 1;
  c1->rc[0][nimd - 1] = c1->ros[0] * sconf->get_rcf(0, nimd - 1);
  extr = c1->rc[0][nimd - 1];
  const double pig = acos(0.0) * 2.0;
  c1->gcs = pig * extr * extr;
  
#ifdef USE_MAGMA
  proc_device = mpidata->rank % device_count;
#else
  proc_device = 0;
#endif

  // the first time the object is used, it will be a "first iteration", to ensure data structures are properly initialised.
  is_first_scale = 1;
  // In the first iteration, if refinement is enabled, determine the number of refinement iterations required to arrive at the target accuracy (if achievable in a reasonable number of iterations)
  refinemode = 2;
  // maxrefiters and accuracygoal should be configurable and preferably set somewhere else
  maxrefiters = 20;
  accuracygoal = 1e-6;
}

InclusionIterationData::InclusionIterationData(const InclusionIterationData& rhs) {
  c1 = new ParticleDescriptorInclusion(reinterpret_cast<ParticleDescriptorInclusion &>(*(rhs.c1)));
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
  vec_zpv = new double[c1->lm * 12];
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double **[12];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[4];
      zpv[zi][zj][0] = vec_zpv + (zi * 12) + (zj * 4);
      zpv[zi][zj][1] = vec_zpv + (zi * 12) + (zj * 4) + 2;
      zpv[zi][zj][0][0] = rhs.zpv[zi][zj][0][0];
      zpv[zi][zj][0][1] = rhs.zpv[zi][zj][0][1];
      zpv[zi][zj][1][0] = rhs.zpv[zi][zj][1][0];
      zpv[zi][zj][1][1] = rhs.zpv[zi][zj][1][1];
    }
  }
  am_vector = new dcomplex[c1->ndm * c1->ndm];
  for (int ai = 0; ai < c1->ndm * c1->ndm; ai++) am_vector[ai] = rhs.am_vector[ai];
  am = new dcomplex*[c1->ndm];
  for (int ai = 0; ai < c1->ndm; ai++) {
    am[ai] = (am_vector + ai * c1->ndm);
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

  nimd = rhs.nimd;
  extr = rhs.extr;
  
  proc_device = rhs.proc_device;
  is_first_scale = rhs.is_first_scale;
  refinemode = rhs.refinemode;
  maxrefiters = rhs.maxrefiters;
  accuracygoal = rhs.accuracygoal;
}

#ifdef MPI_VERSION
InclusionIterationData::InclusionIterationData(const mixMPI *mpidata, const int device_count) {
  c1 = new ParticleDescriptorInclusion(mpidata);
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
  vec_zpv = new double[c1->lm * 12];
  MPI_Bcast(vec_zpv, c1->lm * 12, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  zpv = new double***[c1->lm];
  for (int zi = 0; zi < c1->lm; zi++) {
    zpv[zi] = new double **[12];
    for (int zj = 0; zj < 3; zj++) {
      zpv[zi][zj] = new double*[4];
      zpv[zi][zj][0] = vec_zpv + (zi * 12) + (zj * 4);
      zpv[zi][zj][1] = vec_zpv + (zi * 12) + (zj * 4) + 2;
    }
  }
  am_vector = new dcomplex[c1->ndm * c1->ndm];
  am = new dcomplex*[c1->ndm];
  for (int ai = 0; ai < c1->ndm; ai++) {
    am[ai] = (am_vector + ai * c1->ndm);
    MPI_Bcast(am[ai], c1->ndm, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
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

  MPI_Bcast(&nimd, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&extr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
#ifdef USE_MAGMA
  proc_device = mpidata->rank % device_count;
#else
  proc_device = 0;
#endif
  MPI_Bcast(&refinemode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_first_scale, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&maxrefiters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&accuracygoal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void InclusionIterationData::mpibcast(const mixMPI *mpidata) {
  c1->mpibcast(mpidata);
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
  MPI_Bcast(vec_zpv, c1->lm * 12, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // since MPI expects an int argument for the number of elements to transfer in one go, transfer am one row at a time
  for (int ai = 0; ai < c1->ndm; ai++) {
    MPI_Bcast(am[ai], c1->ndm, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
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
  MPI_Bcast(&nimd, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&extr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&refinemode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_first_scale, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&maxrefiters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&accuracygoal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
#endif

InclusionIterationData::~InclusionIterationData() {
  const int nsph = c1->nsph;
  delete[] am_vector;
  delete[] am;
  for (int zi = 0; zi < c1->lm; zi++) {
    for (int zj = 0; zj < 3; zj++) {
      delete[] zpv[zi][zj];
    }
    delete[] zpv[zi];
  }
  delete[] zpv;
  delete[] vec_zpv;
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

int InclusionIterationData::update_orders(double **rcf, int inner_order, int outer_order) {
  int result = 0;
  int old_lm = c1->lm;
  ((ParticleDescriptorInclusion *)c1)->update_orders(inner_order, outer_order);
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
  instr(rcf, c1);
  thdps(c1->lm, zpv);
  delete[] am;
  delete[] am_vector;
  am_vector = new dcomplex[c1->ndm * c1->ndm]();
  am = new dcomplex*[c1->ndm];
  for (int ai = 0; ai < c1->ndm; ai++) {
    am[ai] = am_vector + ai * c1->ndm;
  }
  return result;
}
// >>> END OF InclusionIterationData CLASS IMPLEMENTATION <<<
