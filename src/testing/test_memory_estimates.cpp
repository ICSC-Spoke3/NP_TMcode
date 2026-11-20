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

/*! \file test_memory_estimates.cpp
 *
 * \brief Check the `get_size()` memory estimates against the real allocations.
 *
 * The program replaces the global `new` and `delete` operators with versions
 * that count the bytes in use, builds the data classes of the calculations
 * from the test input files and compares the bytes that each constructor
 * allocated with the `get_size()` estimate of the class. An estimate passes
 * if it differs from the measured size by at most 1% or 16 KiB, whichever is
 * larger. The program returns 0 if all the estimates pass, 1 otherwise.
 *
 * Usage: `test_memory_estimates [TEST_DATA_FOLDER]`, where the folder
 * defaults to `../../test_data`.
 *
 * \ingroup tests
 */

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <hdf5.h>
#include <new>
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

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

#ifndef INCLUDE_ITERATION_DATA_H_
#include "../include/IterationData.h"
#endif

#ifndef INCLUDE_TFRFME_H_
#include "../include/tfrfme.h"
#endif

using namespace std;

//! \cond
// Heap accounting: each block carries its size in a 16-byte header.
static long live_bytes = 0;
static bool counting = false;

// The replacements are compiled for the host only: in the device pass of an
// OpenMP offload build, clang's device headers define their own operators.
#if !defined(__NVPTX__) && !defined(__AMDGCN__)
void *operator new(size_t size) {
  void *block = malloc(size + 16);
  if (block == NULL) throw std::bad_alloc();
  *(size_t *)block = size;
  if (counting) live_bytes += size;
  return (char *)block + 16;
}
void *operator new[](size_t size) { return operator new(size); }
void operator delete(void *ptr) noexcept {
  if (ptr == NULL) return;
  char *block = (char *)ptr - 16;
  if (counting) live_bytes -= *(size_t *)block;
  free(block);
}
void operator delete[](void *ptr) noexcept { operator delete(ptr); }
void operator delete(void *ptr, size_t) noexcept { operator delete(ptr); }
void operator delete[](void *ptr, size_t) noexcept { operator delete(ptr); }
#endif // host compilation

static int failures = 0;

// Run a constructor call and return the bytes it left allocated.
template <class F> long measure(F construct) {
  long before = live_bytes;
  counting = true;
  construct();
  counting = false;
  return live_bytes - before;
}

static void check(const string &name, long measured, long estimated) {
  long difference = (estimated > measured) ? estimated - measured : measured - estimated;
  long tolerance = measured / 100;
  if (tolerance < 16384) tolerance = 16384;
  bool passed = (difference <= tolerance);
  if (!passed) failures++;
  printf(
    "  %-34s measured %14ld B  estimated %14ld B  %s\n",
    name.c_str(), measured, estimated, passed ? "OK" : "FAILED"
  );
}
//! \endcond

/*! \brief Main program of the memory estimate test.
 *
 * \param argc: `int` Number of command line arguments.
 * \param argv: `char **` Command line arguments (optional test data folder).
 * \return result: `int` 0 if all the estimates pass, 1 otherwise.
 */
int main(int argc, char **argv) {
#ifdef USE_MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  mixMPI *mpidata = new mixMPI(MPI_COMM_WORLD);
#else
  mixMPI *mpidata = new mixMPI();
#endif
  string data = (argc > 1) ? argv[1] : "../../test_data";
  try {
    // Sphere
    {
      ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(data + "/sphere/DEDFB");
      GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(data + "/sphere/DSPH");
      printf("sphere (%s/sphere):\n", data.c_str());
      ParticleDescriptorSphere *pd;
      SphereIterationData *sid, *sid_copy;
      SphereOutputInfo *out, *out_scale;
      check("ParticleDescriptorSphere", measure([&] { pd = new ParticleDescriptorSphere(gconf, sconf); }),
	    ParticleDescriptorSphere::get_size(gconf, sconf));
      check("SphereIterationData", measure([&] { sid = new SphereIterationData(gconf, sconf, mpidata, 0); }),
	    SphereIterationData::get_size(gconf, sconf));
      check("SphereIterationData (copy)", measure([&] { sid_copy = new SphereIterationData(*sid); }),
	    SphereIterationData::get_size(gconf, sconf));
      check("SphereOutputInfo", measure([&] { out = new SphereOutputInfo(sconf, gconf, mpidata); }),
	    SphereOutputInfo::get_size(sconf, gconf));
      check("SphereOutputInfo (one scale)", measure([&] { out_scale = new SphereOutputInfo(sconf, gconf, mpidata, 1, 1); }),
	    SphereOutputInfo::get_size(sconf, gconf, 1, 1));
      delete pd; delete sid; delete sid_copy; delete out; delete out_scale;
      delete sconf; delete gconf;
    }
    // Cluster
    {
      ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(data + "/cluster/DEDFB");
      GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(data + "/cluster/DCLU");
      printf("cluster (%s/cluster):\n", data.c_str());
      ParticleDescriptorCluster *pd;
      ClusterIterationData *cid, *cid_copy;
      ClusterOutputInfo *out, *out_scale;
      check("ParticleDescriptorCluster", measure([&] { pd = new ParticleDescriptorCluster(gconf, sconf); }),
	    ParticleDescriptorCluster::get_size(gconf, sconf));
      check("ClusterIterationData", measure([&] { cid = new ClusterIterationData(gconf, sconf, mpidata, 0); }),
	    ClusterIterationData::get_size(gconf, sconf));
      check("ClusterIterationData (copy)", measure([&] { cid_copy = new ClusterIterationData(*cid); }),
	    ClusterIterationData::get_size(gconf, sconf));
      check("ClusterOutputInfo", measure([&] { out = new ClusterOutputInfo(sconf, gconf, mpidata); }),
	    ClusterOutputInfo::get_size(sconf, gconf));
      check("ClusterOutputInfo (one scale)", measure([&] { out_scale = new ClusterOutputInfo(sconf, gconf, mpidata, 1, 1); }),
	    ClusterOutputInfo::get_size(sconf, gconf, 1, 1));
      delete pd; delete cid; delete cid_copy; delete out; delete out_scale;
      delete sconf; delete gconf;
    }
    // Inclusion
    {
      ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(data + "/inclusion/DEDFB");
      GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(data + "/inclusion/DINCLU");
      printf("inclusion (%s/inclusion):\n", data.c_str());
      ParticleDescriptorInclusion *pd;
      InclusionIterationData *iid, *iid_copy;
      InclusionOutputInfo *out, *out_scale;
      check("ParticleDescriptorInclusion", measure([&] { pd = new ParticleDescriptorInclusion(gconf, sconf); }),
	    ParticleDescriptorInclusion::get_size(gconf, sconf));
      check("InclusionIterationData", measure([&] { iid = new InclusionIterationData(gconf, sconf, mpidata, 0); }),
	    InclusionIterationData::get_size(gconf, sconf));
      check("InclusionIterationData (copy)", measure([&] { iid_copy = new InclusionIterationData(*iid); }),
	    InclusionIterationData::get_size(gconf, sconf));
      check("InclusionOutputInfo", measure([&] { out = new InclusionOutputInfo(sconf, gconf, mpidata); }),
	    InclusionOutputInfo::get_size(sconf, gconf));
      check("InclusionOutputInfo (one scale)", measure([&] { out_scale = new InclusionOutputInfo(sconf, gconf, mpidata, 1, 1); }),
	    InclusionOutputInfo::get_size(sconf, gconf, 1, 1));
      delete pd; delete iid; delete iid_copy; delete out; delete out_scale;
      delete sconf; delete gconf;
    }
    // Trapping (sizes of the performance test case, with a reduced NKV)
    {
      const int lm = 20, nkv = 21, nv = 21;
      printf("trapping (LM = %d, NKV = %d, NXV = NYV = NZV = %d):\n", lm, nkv, nv);
      Swap1 *swap1;
      Swap2 *swap2;
      TFRFME *tfrfme;
      TrappingOutputInfo *out;
      check("Swap1", measure([&] { swap1 = new Swap1(lm, nkv); }), Swap1::get_size(lm, nkv));
      check("Swap2", measure([&] { swap2 = new Swap2(nkv); }), Swap2::get_size(nkv));
      check("TFRFME", measure([&] { tfrfme = new TFRFME(0, lm, nkv, nv, nv, nv); }),
	    TFRFME::get_size(lm, nkv, nv, nv, nv));
      check("TrappingOutputInfo", measure([&] { out = new TrappingOutputInfo(0, 1, 0, nv, nv, nv); }),
	    TrappingOutputInfo::get_size(0, 1, nv, nv, nv));
      delete swap1; delete swap2; delete tfrfme; delete out;
    }
  } catch (const exception &ex) {
    printf("[ERROR] %s\n", ex.what());
    failures++;
  }
  printf("%s: %d estimate(s) outside the tolerance.\n", (failures == 0) ? "SUCCESS" : "FAILURE", failures);
  delete mpidata;
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return (failures == 0) ? 0 : 1;
}
