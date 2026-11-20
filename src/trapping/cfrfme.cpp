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

/*! \file cfrfme.cpp
 *
 * C++ implementation of FRFME functions.
 * \ingroup prog_trapping
 */
#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <regex>
#include <string>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_PARSERS_H_
#include "../include/Parsers.h"
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

#ifndef INCLUDE_TFRFME_H_
#include "../include/tfrfme.h"
#endif

#ifndef INCLUDE_TRA_SUBS_H_
#include "../include/tra_subs.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif

#ifndef INCLUDE_RESOURCES_H_
#include "../include/resources.h"
#endif

#ifdef USE_NVTX
#include <nvtx3/nvToolsExt.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_TARGET_OFFLOAD
#include <cstdlib>

#ifndef FRFME_ORDER_TILE
/*! \brief Number of orders summed together by each GPU thread in offload_loop().
 *
 * Each GPU thread keeps the sums of this many consecutive orders in
 * registers while it loops over the wave vectors, and writes them to the
 * vector of weight sums once at the end. Larger values compute the phase
 * factors fewer times, but use more registers per thread. The default can
 * be changed at build time with `-DFRFME_ORDER_TILE=<n>`.
 */
#define FRFME_ORDER_TILE 16
#endif

/*! \brief Specialized function to perform GPU-offloaded trapping loop.
 *
 * The offload of GPU operations through interface layers, such as OpenMP,
 * can become implementation dependent, especially when attempting to target
 * generic hardware. In this case, not all the compilers are able to obtain
 * a correct solution, especially when loop optimization is requested. In
 * order to preserve the work-flow integrity when compiling in high optimization,
 * therefore, it is useful to encapsulate the GPU oflloaded loop in a separate
 * function, which takes care of bracketing the target code region allowing
 * for the loop execution on a generic target device.
 *
 * \param vec_wsum: `dcomplex *` Vector of weight sums to be computed.
 * \param size_vec_wsum: `int` Size of the vector of weight sums to be computed.
 * \param global_vec_w: `dcomplex *` Global space for weight vectors visible to
 * all threads.
 * \param size_global_vec_w: `const int` Size of the global weight vector space.
 * \param vec_tt1_wk: `dcomplex *` Vector of swap object weights.
 * \param size_vec_tt1_wk: `const int` Size of the swap object weights vector.
 * \param vkv: `double *` Vector of beam fields wave numbers (with size = NLMMT).
 * \param _xv: `double *` Vector of x-coordinate calculation spacing.
 * \param nxv: `int` Number of x-coordinates used in the calculation.
 * \param _yv: `double *` Vector of y-coordinate calculation spacing.
 * \param nyv: `int` Number of y-coordinates used in the calculation.
 * \param _zv: `double *` Vector of z-coordinate calculation spacing.
 * \param nzv: `int` Number of z-coordinates used in the calculation.
 * \param vec_vkzm: `double *` Vectorized VKZM matrix (size = NLMMT x NLMMT)
 * \param jlmf: `int` First order of the calculation.
 * \param jlml: `int` Last order of the calculation.
 * \param nkv: `int` Number of beam vector wave-numbers.
 * \param nlmmt: `int` NLMMT = 2 x LM x (LM + 2), with LM as maximum field
 * expansion order.
 * \param delks: `double` Square of the wave-number grid spacing.
 * \param frsh: `double` Instrumental offset along beam z-axis.
 */
void offload_loop(
  dcomplex *vec_wsum, int size_vec_wsum, dcomplex *global_vec_w, const int size_global_vec_w,
  const dcomplex *vec_tt1_wk, const int size_vec_tt1_wk, double *vkv, double *_xv, const int nxv,
  double *_yv, const int nyv, double *_zv, const int nzv, double *vec_vkzm, const int jlmf, const int jlml,
  const int nkv, const int nlmmt, double delks, double frsh
);
#endif

#ifndef USE_TARGET_OFFLOAD
#ifndef FRFME_CPU_ORDER_TILE
/*! \brief Number of orders summed together by each thread in cpu_loop().
 *
 * The phase factors of a block of grid points are computed once for each
 * tile of this many orders. Can be changed at build time with
 * `-DFRFME_CPU_ORDER_TILE=<n>`.
 */
#define FRFME_CPU_ORDER_TILE 32
#endif

#ifndef FRFME_CPU_POINT_BLOCK
/*! \brief Number of grid points handled together by each thread in cpu_loop().
 *
 * Each coefficient of the beam expansion is read once for a block of this
 * many grid points. Can be changed at build time with
 * `-DFRFME_CPU_POINT_BLOCK=<n>`.
 */
#define FRFME_CPU_POINT_BLOCK 16
#endif

/*! \brief Compute the weight sums of the trapping loop on the CPU.
 *
 * For every order \f$j\f$ from `jlmf` to `jlml` and every grid point, this
 * function computes the sum over the beam wave vectors of the beam
 * expansion coefficient of order \f$j\f$ times the phase factor of the
 * wave vector at the grid point, and stores it, multiplied by `delks`, in
 * `vec_wsum`. The rows of the other orders are left untouched.
 *
 * The work is split in tiles of FRFME_CPU_ORDER_TILE orders and blocks of
 * FRFME_CPU_POINT_BLOCK grid points. Each thread computes the phase factors
 * of its block of grid points once per tile of orders, and keeps the sums
 * of the tile in a small local array, so that the phase factors are not
 * recomputed for every order and each coefficient is read once per block
 * of grid points. The sums over the wave vectors run in the same order as
 * in offload_loop().
 *
 * \param vec_wsum: `dcomplex *` Vector of weight sums (size = NLMMT x NXV x NYV x NZV).
 * \param vec_tt1_wk: `const dcomplex *` Vector of beam expansion coefficients
 * (size = NKV x NKV x NLMMT).
 * \param vkv: `double *` Vector of wave-number values.
 * \param _xv: `double *` Vector of x-coordinates of the grid.
 * \param nxv: `int` Number of x-coordinates.
 * \param _yv: `double *` Vector of y-coordinates of the grid.
 * \param nyv: `int` Number of y-coordinates.
 * \param _zv: `double *` Vector of z-coordinates of the grid.
 * \param nzv: `int` Number of z-coordinates.
 * \param vec_vkzm: `double *` Vectorized VKZM matrix (size = NKV x NKV).
 * \param jlmf: `int` First order of the calculation.
 * \param jlml: `int` Last order of the calculation.
 * \param nkv: `int` Number of beam vector wave-numbers.
 * \param nlmmt: `int` NLMMT = 2 x LM x (LM + 2), with LM as maximum field
 * expansion order.
 * \param delks: `double` Square of the wave-number grid spacing.
 * \param frsh: `double` Instrumental offset along beam z-axis.
 */
void cpu_loop(
  dcomplex *vec_wsum, const dcomplex *vec_tt1_wk, double *vkv, double *_xv,
  const int nxv, double *_yv, const int nyv, double *_zv, const int nzv,
  double *vec_vkzm, const int jlmf, const int jlml, const int nkv,
  const int nlmmt, double delks, double frsh
);
#endif // USE_TARGET_OFFLOAD

using namespace std;

//! \cond
// Index of the line of the input file being parsed, for error messages.
static int frfme_input_line = 0;
//! \endcond

/*! \brief C++ implementation of FRFME, without the input error handling.
 *
 *  \param data_file: `string` Name of the input data file.
 *  \param output_path: `string` Directory to write the output files in.
 */
static void frfme_unguarded(string data_file, string output_path);

/*! \brief C++ implementation of FRFME
 *
 *  Malformed values in the input file are reported as a `ParseError` that
 *  names the file and the line.
 *
 *  \param data_file: `string` Name of the input data file.
 *  \param output_path: `string` Directory to write the output files in.
 */
void frfme(string data_file, string output_path) {
  frfme_input_line = 0;
  parse_guarded(data_file, frfme_input_line, [&]() { frfme_unguarded(data_file, output_path); });
}

static void frfme_unguarded(string data_file, string output_path) {
#ifdef USE_NVTX
  nvtxRangePush("Running frfme()");
#endif
  chrono::time_point<chrono::high_resolution_clock> t_start = chrono::high_resolution_clock::now();
  chrono::time_point<chrono::high_resolution_clock> t_end;
  chrono::duration<double> elapsed, frfme_duration;
  char buffer[256];
  string message = "INIT";
  Logger logger(LOG_INFO);
  string tfrfme_name = output_path + "/c_TFRFME.hd5";
  TFRFME *tfrfme = NULL;
  Swap1 *tt1 = NULL;
  Swap2 *tt2 = NULL;
  char namef[7];
  char more;
  dcomplex *wk = NULL;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  int line_count = 0, last_read_line = 0;
  regex re = regex("-?[0-9]+");
  string *file_lines = load_file(data_file, &line_count);
  smatch m;
  string str_target = checked_line(file_lines, line_count, last_read_line++, data_file, frfme_input_line);
  regex_search(str_target, m, re);
  int jlmf = stoi(m.str());
  str_target = m.suffix().str();
  regex_search(str_target, m, re);
  int jlml = stoi(m.str());
  // Input parameters
  int lmode, lm;
  double vk = 0.0, exri = 0.0, an = 0.0, ff = 0.0, tra = 0.0;
  double exdc = 0.0, wp = 0.0, xip = 0.0, xi = 0.0;
  int idfc = 0, nxi = 0;
  double apfafa = 0.0, pmf = 0.0, spd = 0.0, rir = 0.0, ftcn = 0.0, fshmx = 0.0;
  double vxyzmx = 0.0, delxyz = 0.0, vknmx = 0.0, delk = 0.0, delks = 0.0, vkm = 0.0;
  double frsh = 0.0, exril = 1.0, spdfr = 1.0, exdcl = 1.0, wlenfr= 1.0, wn = 0.0;
  int nksh = 0, nrsh = 0, nxsh = 0, nysh = 0, nzsh = 0;
  // Calculation parameters, computed from input data
  int nlmmt = 0, nrvc = 0;
  int nkshpo = 0, nks = 0, nkv = 0;
  int nxshpo = 0, nxs = 0, nxv = 0;
  int nyshpo = 0, nys = 0, nyv = 0;
  int nzshpo = 0, nzs = 0, nzv = 0;
  double *vkv = NULL;
  double *_xv = NULL, *_yv = NULL, *_zv = NULL;
  // Vector size variables
  int wsum_size;
  // End of vector size variables
  if (jlmf != 1) {
#ifdef USE_NVTX
    nvtxRangePush("frfme() with jlmf != 1");
#endif
    if (tfrfme == NULL) tfrfme = TFRFME::from_binary(tfrfme_name, "HDF5");
    if (tfrfme != NULL) {
      lmode = tfrfme->lmode;
      lm = tfrfme->lm;
      nkv = tfrfme->nkv;
      nks = nkv - 1;
      nxv = tfrfme->nxv;
      nyv = tfrfme->nyv;
      nzv = tfrfme->nzv;
      _xv = tfrfme->get_x();
      _yv = tfrfme->get_y();
      _zv = tfrfme->get_z();
      vk = tfrfme->vk;
      exri = tfrfme->exri;
      an = tfrfme->an;
      ff = tfrfme->ff;
      tra = tfrfme->tra;
      spd = tfrfme->spd;
      frsh = tfrfme->frsh;
      exril = tfrfme->exril;
      string tempname2 = output_path + "/c_TEMPTAPE2.hd5";
      if (tt2 == NULL) tt2 = Swap2::from_binary(tempname2, "HDF5");
      if (tt2 != NULL) {
	vkv = tt2->get_vector();
	apfafa = tt2->apfafa;
	pmf = tt2->pmf;
	spd = tt2->spd;
	rir = tt2->rir;
	ftcn = tt2->ftcn;
	fshmx = tt2->fshmx;
	vxyzmx = tt2->vxyzmx;
	delxyz = tt2->delxyz;
	vknmx = tt2->vknmx;
	delk = tt2->delk;
	delks = tt2->delks;
	nlmmt = tt2->nlmmt;
	nrvc = tt2->nrvc;
      } else {
	message = "could not open TEMPTAPE2 file.\n";
	logger.err(message);
      }
      string tempname1 = output_path + "/c_TEMPTAPE1.hd5";
      if (tt1 == NULL) tt1 = Swap1::from_binary(tempname1, "HDF5");
      if (tt1 == NULL) {
	message = "could not open TEMPTAPE1 file.\n";
	logger.err(message);
      }
    } else {
      message = "could not open TFRFME file.\n";
      logger.err(message);
    }
#ifdef USE_NVTX
    nvtxRangePop();
#endif
  } else { // label 16; jlmf = 1
#ifdef USE_NVTX
    nvtxRangePush("frfme() with jlmf == 1");
#endif
#ifdef USE_NVTX
    nvtxRangePush("Setup operations");
#endif
    str_target = checked_line(file_lines, line_count, last_read_line++, data_file, frfme_input_line);
    for (int cli = 0; cli < 7; cli++) {
      regex_search(str_target, m, re);
      if (cli == 0) lmode = stoi(m.str());
      else if (cli == 1) lm = stoi(m.str());
      else if (cli == 2) nksh = stoi(m.str());
      else if (cli == 3) nrsh = stoi(m.str());
      else if (cli == 4) nxsh = stoi(m.str());
      else if (cli == 5) nysh = stoi(m.str());
      else if (cli == 6) nzsh = stoi(m.str());
      str_target = m.suffix().str();
    }
    re = regex("-?[0-9]\\.[0-9]+([dDeE][-+]?[0-9]+)?");
    regex_search(str_target, m, re);
    wlenfr = stod(m.str());
    str_target = checked_line(file_lines, line_count, last_read_line++, data_file, frfme_input_line);
    for (int cli = 0; cli < 3; cli++) {
      regex_search(str_target, m, re);
      if (cli == 0) an = stod(m.str());
      else if (cli == 1) ff = stod(m.str());
      else if (cli == 2) tra = stod(m.str());
      str_target = m.suffix().str();
    }
    str_target = checked_line(file_lines, line_count, last_read_line++, data_file, frfme_input_line);
    for (int cli = 0; cli < 3; cli++) {
      regex_search(str_target, m, re);
      if (cli == 0) spd = stod(m.str());
      else if (cli == 1) spdfr = stod(m.str());
      else if (cli == 2) exdcl = stod(m.str());
      str_target = m.suffix().str();
    }
    str_target = checked_line(file_lines, line_count, last_read_line++, data_file, frfme_input_line);
    re = regex("[eEmM]");
#ifdef USE_NVTX
    nvtxRangePop();
#endif
    if (regex_search(str_target, m, re)) {
      more = m.str().at(0);
      if (more == 'm' || more == 'M') {
	more = 'M';
	sprintf(namef, "c_TMDF");
      }
      else if (more == 'e' || more == 'E') {
	more = 'E';
	sprintf(namef, "c_TEDF");
      } else {
	throw(UnrecognizedConfigurationException("[ERROR] only 'm', 'M', 'e', or 'E' accepted as modes!\n"));
      }
    }
    str_target = m.suffix().str();
    re = regex("[0-9]+");
    regex_search(str_target, m, re);
    int ixi = stoi(m.str());
    string tedf_name = output_path + "/" + namef + ".hd5";
    // Check for run-time options
    bool skip_frfme = false;
    double host_ram_gb = 0.0;
    short memory_policy = MEMORY_POLICY_STOP;
    last_read_line++; // skip the end-of-data code
    while (last_read_line < line_count) {
      bool is_parsed = false;
      str_target = checked_line(file_lines, line_count, last_read_line++, data_file, frfme_input_line);
      if (str_target.size() > 14 && str_target.substr(0, 14).compare("MEMORY_POLICY=") == 0) {
	memory_policy = parse_memory_policy(str_target.substr(14));
	is_parsed = true;
      }
      if (str_target.size() > 12 && str_target.substr(0, 12).compare("HOST_RAM_GB=") == 0) {
	host_ram_gb = stod(str_target.substr(12));
	is_parsed = true;
      }
      if (str_target.size() > 12) {
	if (str_target.substr(0, 12).compare("PRECOMPUTED=") == 0) {
	  string precomp_name = output_path + "/" + str_target.substr(12, str_target.size());
	  // TODO: check that the file exists and it is suitable for the calculation.
	  message = "using precomputed file " + precomp_name + "\n";
	  logger.log(message);
	  skip_frfme = true;
	  is_parsed = true;
	}
      }
      if (!is_parsed) {
	if (str_target.size() > 0) {
	  if (str_target.substr(0, 1).compare("#") != 0) {
	    // ERROR: found an unrecognized option that is not marked as comment
	    throw(UnrecognizedConfigurationException("[ERROR] unrecognized option \"" + str_target + "\"!\n"));
	  }
	}
      }
    }
    if (skip_frfme) return;
    check_input_file(tedf_name, "scatterer data (written by a previous np_sphere or np_cluster run)");
    ScattererConfiguration *tedf = ScattererConfiguration::from_binary(tedf_name, "HDF5");
#ifdef USE_NVTX
    nvtxRangePush("TEDF data import");
#endif
    if (tedf == NULL) {
      throw(InputError("could not open the scatterer data file \"" + tedf_name + "\""));
    }
    int iduml, idum;
    iduml = tedf->number_of_spheres;
    idum = tedf->get_iog(iduml - 1);
    exdc = tedf->exdc;
    wp = tedf->wp;
    xip = tedf->xip;
    idfc = tedf->idfc;
    nxi = tedf->number_of_scales;
    if (idfc >= 0) {
      if (ixi <= nxi) {
	xi = tedf->get_scale(ixi - 1);
      } else { // label 96
	delete tedf;
	// label 98
	string output_name = output_path + "/c_OFRFME";
	FILE *output = fopen(output_name.c_str(), "w");
	fprintf(output, "  WRONG INPUT TAPE\n");
	fclose(output);
      }
    } else { // label 18
      xi = xip;
    }
    // label 20
#ifdef USE_NVTX
    nvtxRangePop();
#endif
    delete tedf;
    wn = wp / 3.0e8;
    vk = xi * wn;
    exri = sqrt(exdc);
    frsh = 0.0;
    exril = 0.0;
    fshmx = 0.0;
    apfafa = exri / (an * ff);
    if (lmode != 0) pmf = 2.0 * apfafa;
    if (spd > 0.0) {
      exril = sqrt(exdcl);
      rir = exri / exril;
      ftcn = 2.0 / (1.0 + rir);
      frsh = -spd * spdfr;
      double sthmx = an / exri;
      double sthlmx = sthmx * rir;
      double uy = 1.0;
      fshmx = spd * (rir * (sqrt(uy - sthmx * sthmx) / sqrt(uy - sthlmx * sthlmx)) - uy);
    }
    // label 22
#ifdef USE_NVTX
    nvtxRangePush("Memory data loading");
#endif
    nlmmt = lm * (lm + 2) * 2;
    nks = nksh * 2;
    nkv = nks + 1;
    // Array initialization
    long swap1_size, swap2_size, tfrfme_size;
    double size_mb;
    message = "calculating memory requirements...\n";
    logger.log(message);
    swap1_size = Swap1::get_size(lm, nkv);
    size_mb = 1.0 * swap1_size / 1024.0 / 1024.0;
    sprintf(buffer, "Swap 1: %.2lg MB\n", size_mb);
    message = string(buffer);
    logger.log(message);
    swap2_size = Swap2::get_size(nkv);
    size_mb = 1.0 * swap2_size / 1024.0 / 1024.0;
    sprintf(buffer, "Swap 2: %.2lg MB\n", size_mb);
    message = string(buffer);
    logger.log(message);
    tt2 = new Swap2(nkv);
    vkv = tt2->get_vector();
    // End of array initialization
    vkm = vk * exri;
    vknmx = vk * an;
    delk = vknmx / nksh;
    delks = delk / vkm;
    delks = delks * delks;
    vxyzmx = acos(0.0) * 4.0 / vkm * wlenfr;
    delxyz = vxyzmx / nrsh;
    nxs = nxsh * 2;
    nxv = nxs + 1;
    nxshpo = nxsh + 1;
    nys = nysh * 2;
    nyv = nys + 1;
    nyshpo = nysh + 1;
    nzs = nzsh * 2;
    nzv = nzs + 1;
    nzshpo = nzsh + 1;
    tfrfme_size = TFRFME::get_size(lm, nkv, nxv, nyv, nzv);
    size_mb = 1.0 * tfrfme_size / 1024.0 / 1024.0;
    sprintf(buffer, "TFRFME: %.2lg MB\n", size_mb);
    message = string(buffer);
    logger.log(message);
    // Coefficients rearranged for the loop over grid points: by tiles of
    // orders in cpu_loop(), or copied for the GPU in offload_loop().
    const long loop_nkvs = (long)nkv * nkv;
    const long loop_orders = jlml - jlmf + 1;
#ifdef USE_TARGET_OFFLOAD
    long loop_size = (long)sizeof(dcomplex) * loop_nkvs * loop_orders;
    loop_size = ((loop_size + 63) / 64) * 64;
#else
    const long loop_tiles = (loop_orders + FRFME_CPU_ORDER_TILE - 1) / FRFME_CPU_ORDER_TILE;
    long loop_size = (long)sizeof(dcomplex) * loop_tiles * FRFME_CPU_ORDER_TILE * loop_nkvs;
    loop_size += 4 * (long)sizeof(double) * loop_nkvs; // wave vector components and weights
#endif // USE_TARGET_OFFLOAD
    size_mb = 1.0 * loop_size / 1024.0 / 1024.0;
    sprintf(buffer, "Loop coefficients: %.2lg MB\n", size_mb);
    message = string(buffer);
    logger.log(message);
    size_mb = 1.0 * (swap1_size + swap2_size + tfrfme_size + loop_size) / 1024.0 / 1024.0;
    sprintf(buffer, "TOTAL: %.2lg MB\n", size_mb);
    message = string(buffer);
    logger.log(message);
    // Swap 2 is already allocated: check the rest before allocating it.
    check_host_memory(
      swap1_size + tfrfme_size + loop_size, "the FRFME data", host_ram_gb, memory_policy, &logger
    );
    sprintf(buffer, "nxv = %d, nyv = %d, nzv = %d, nrvc = %d\n", nxv, nyv, nzv, (nxv * nyv * nzv));
    message = string(buffer);
    logger.log(message);
    sprintf(buffer, "nlmmt = %d, nkv = %d\n", nlmmt, nkv);
    message = string(buffer);
    logger.log(message);
    
    tfrfme = new TFRFME(lmode, lm, nkv, nxv, nyv, nzv);
    if (tfrfme == NULL)
      throw(ObjectAllocationException("[ERROR] could not initialize TFRFME structure!\n"));
    _xv = tfrfme->get_x();
    _yv = tfrfme->get_y();
    _zv = tfrfme->get_z();
#ifdef USE_NVTX
    nvtxRangePop();
#endif
#ifdef USE_NVTX
    nvtxRangePush("Looped vector initialization");
#endif
    for (int i24 = nxshpo; i24 <= nxs; i24++) {
      _xv[i24] = _xv[i24 - 1] + delxyz;
      _xv[nxv - i24 - 1] = -_xv[i24];
    } // i24 loop
    for (int i25 = nyshpo; i25 <= nys; i25++) {
      _yv[i25] = _yv[i25 - 1] + delxyz;
      _yv[nyv - i25 - 1] = -_yv[i25];
    } // i25 loop
    for (int i27 = nzshpo; i27 <= nzs; i27++) {
      _zv[i27] = _zv[i27 - 1] + delxyz;
      _zv[nzv - i27 - 1] = -_zv[i27];
    } // i27 loop
    nrvc = nxv * nyv * nzv;
    nkshpo = nksh + 1;
    for (int i28 = nkshpo; i28 <= nks; i28++) {
      vkv[i28] = vkv[i28 - 1] + delk;
      vkv[nkv - i28 - 1] = -vkv[i28];
    } // i28 loop
#ifdef USE_NVTX
    nvtxRangePop();
#endif
#ifdef USE_NVTX
    nvtxRangePush("TFRFME initialization");
#endif
    tfrfme->set_param("vk", vk);
    tfrfme->set_param("exri", exri);
    tfrfme->set_param("an", an);
    tfrfme->set_param("ff", ff);
    tfrfme->set_param("tra", tra);
    tfrfme->set_param("spd", spd);
    tfrfme->set_param("frsh", frsh);
    tfrfme->set_param("exril", exril);
    //fstream temptape1, temptape2;
    string temp_name1 = output_path + "/c_TEMPTAPE1.hd5";
    string temp_name2 = output_path + "/c_TEMPTAPE2.hd5";
    //temptape1.open(temp_name1.c_str(), ios::out | ios::binary);
    tt1 = new Swap1(lm, nkv);
    if (wk != NULL) delete[] wk;
    wk = frfmer(nkv, vkm, vknmx, apfafa, tra, spd, rir, ftcn, lm, lmode, pmf, tt1, tt2);
    tt1->write_binary(temp_name1, "HDF5");
    //delete tt1;
    tt2->set_param("apfafa", apfafa);
    tt2->set_param("pmf", pmf);
    tt2->set_param("spd", spd);
    tt2->set_param("rir", rir);
    tt2->set_param("ftcn", ftcn);
    tt2->set_param("fshmx", fshmx);
    tt2->set_param("vxyzmx", vxyzmx);
    tt2->set_param("delxyz", delxyz);
    tt2->set_param("vknmx", vknmx);
    tt2->set_param("delk", delk);
    tt2->set_param("delks", delks);
    tt2->set_param("nlmmt", 1.0 * nlmmt);
    tt2->set_param("nrvc", 1.0 * nrvc);
    tt2->write_binary(temp_name2, "HDF5");
#ifdef USE_NVTX
    nvtxRangePop();
#endif
  } // jlmf = 1 configuration mode
  // label 45
#ifdef USE_NVTX
  nvtxRangePush("j80 loop");
#endif
  const dcomplex *vec_tt1_wk = tt1->wk;
  dcomplex *vec_wsum = tfrfme->vec_wsum;
  double *vec_vkzm = tt2->vec_vkzm;
#ifdef USE_TARGET_OFFLOAD
  const int nkvs = nkv * nkv;
  int size_vec_wsum = nlmmt * nrvc;
  int size_global_vec_w = nkvs * (jlml - jlmf + 1);
  int size_vec_tt1_wk = nkvs * nlmmt;
  dcomplex *global_vec_w = (dcomplex *)zeroed_aligned_alloc(64, size_global_vec_w * sizeof(dcomplex));
#endif // USE_TARGET_OFFLOAD
  message = "looping over " + to_string(jlml - jlmf + 1) + " J iterations.\n";
  logger.log(message);
#ifdef USE_TARGET_OFFLOAD
  t_end = chrono::high_resolution_clock::now();
  elapsed = t_start - t_end;
  frfme_duration = elapsed;
  t_start = chrono::high_resolution_clock::now();
  message = "computing loop.\n";
  logger.log(message);
#ifdef USE_NVTX
  nvtxRangePush("Offloaded loop");
#endif
  offload_loop(
    vec_wsum, size_vec_wsum, global_vec_w, size_global_vec_w, vec_tt1_wk,
    size_vec_tt1_wk, vkv, _xv, nxv, _yv, nyv, _zv, nzv, vec_vkzm, jlmf, jlml,
    nkv, nlmmt, delks, frsh
  );
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  t_end = chrono::high_resolution_clock::now();
  elapsed = t_end - t_start;
  frfme_duration += elapsed;
  sprintf(buffer, "loop calculation took %lfs.\n", elapsed.count());
  message = string(buffer);
  logger.log(message);
  free(global_vec_w);
#else
  // This code block is compiled if USE_TARGET_OFFLOAD is not defined
  t_start = chrono::high_resolution_clock::now();
  message = "computing loop.\n";
  logger.log(message);
  cpu_loop(
    vec_wsum, vec_tt1_wk, vkv, _xv, nxv, _yv, nyv, _zv, nzv, vec_vkzm,
    jlmf, jlml, nkv, nlmmt, delks, frsh
  );
  t_end = chrono::high_resolution_clock::now();
  elapsed = t_end - t_start;
  sprintf(buffer, "loop calculation took %lfs.\n", elapsed.count());
  message = string(buffer);
  logger.log(message);
#endif // USE_TARGET_OFFLOAD
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  // label 88
#ifdef USE_NVTX
  nvtxRangePush("Closing operations");
#endif
  tfrfme->write_binary(tfrfme_name, "HDF5");
  string output_name = output_path + "/c_OFRFME";
  FILE *output = fopen(output_name.c_str(), "w");
  fprintf(output, " IF JLML < NLMMT, PRESERVE TEMPTAPE1, TEMPTAPE2, AND TFRFME,\n");
  fprintf(output, " AND RESTART LM RUN WITH JLMF = JLML+1\n");
  if (spd > 0.0) fprintf(output, "  FSHMX =%15.7lE\n", fshmx);
  fprintf(output, "  FRSH =%15.7lE\n", frsh);
  fclose(output);
#ifdef USE_NVTX
  nvtxRangePop();
  nvtxRangePop();
#endif
  // label 45
#ifdef USE_NVTX
  nvtxRangePush("frfme() memory clean");
#endif
  if (tfrfme != NULL) delete tfrfme;
  delete[] file_lines;
  if (tt2 != NULL) delete tt2;
  if (wk != NULL) delete[] wk;
  if (tt1 != NULL) delete tt1;
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_TARGET_OFFLOAD
  elapsed = chrono::high_resolution_clock::now() - t_end;
  frfme_duration += elapsed;
#else
  elapsed = chrono::high_resolution_clock::now() - t_start;
  frfme_duration = elapsed;
#endif
  message = "FRFME took " + to_string(frfme_duration.count()) + "s.\n";
  logger.log(message);
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}
 // STALE CODE BLOCK
 // else { // Should never happen.
 // 	  message = "ERROR: could not open TFRFME file for output.\n";
 // 	  logger.err(message);
 // 	}
 //      } else {
 // 	message = "ERROR: could not open TEDF file.\n";
 // 	logger.err(message);
 //      }
 //    } else { // label 98
 //      string output_name = output_path + "/c_OFRFME";
 //      FILE *output = fopen(output_name.c_str(), "w");
 //      fprintf(output, "  WRONG INPUT TAPE\n");
 //      fclose(output);
 //    }
 // #ifdef USE_NVTX
 //     nvtxRangePop();
 // #endif
 //   }

#ifdef USE_TARGET_OFFLOAD
void offload_loop(
  dcomplex *vec_wsum, const int size_vec_wsum, dcomplex *global_vec_w, const int size_global_vec_w,
  const dcomplex *vec_tt1_wk, const int size_vec_tt1_wk, double *vkv, double *_xv, const int nxv,
  double *_yv, const int nyv, double *_zv, const int nzv, double *vec_vkzm, const int jlmf, const int jlml,
  const int nkv, const int nlmmt, double delks, double frsh
) {
  long long nvtot = nxv * nyv * nzv;
  long long nkvs = nkv * nkv;
  long long nkvmo = nkv - 1; 
  long long nkvvmo = nkvmo * nkv;
  long long nvxy = nxv * nyv;
  dcomplex cc0 = 0.0 + I * 0.0;

 // Inizializza global_vec_w e vec_wsum sulla CPU in parallelo
#pragma omp parallel for collapse(2)
  for (int j80 = jlmf - 1; j80 < jlml; j80++) {
    for (long long jxy50 = 0; jxy50 < nkvs; jxy50++) {
      int j80_index = j80 - jlmf + 1;
      dcomplex *vec_w = global_vec_w + nkvs * j80_index;
      long long wk_index = nlmmt * jxy50;
      dcomplex wk_value = vec_tt1_wk[wk_index + j80];
      long long jy50 = jxy50 / nkv;
      long long jx50 = jxy50 % nkv;
      vec_w[(nkv * jx50) + jy50] = wk_value;
    } // jxy50 loop
  }
  
  // Clear only the rows of the orders computed in this run: when resuming
  // (jlmf > 1), the rows of the earlier orders hold results read from TFRFME.
#pragma omp parallel for
  for (long long i = (jlmf - 1) * nvtot; i < jlml * nvtot; i++) {
    vec_wsum[i] = cc0;
  }

  constexpr int order_tile = FRFME_ORDER_TILE;
  const int num_orders = jlml - jlmf + 1;
  const int num_tiles = (num_orders + order_tile - 1) / order_tile;
  const double delks_value = delks;

#pragma omp target data map(tofrom: vec_wsum[0:size_vec_wsum]) \
  map(to: global_vec_w[0:size_global_vec_w]) \
  map(to: _xv[0:nxv], _yv[0:nyv], _zv[0:nzv]) \
  map(to: vkv[0:nkv], vec_vkzm[0:nkvs])
  {
    // Each GPU thread handles one grid point (ixyz) and one tile of up to
    // FRFME_ORDER_TILE consecutive orders. It sums over the wave vectors
    // (jy60x55) in registers and writes each element of vec_wsum once, at
    // the end. No two threads write the same element, so the accumulation
    // needs neither atomics nor a reduction, and the result does not depend
    // on the number of threads. Writing the sums once, instead of reading and
    // writing vec_wsum in global memory for every term, removes the memory
    // traffic that limited the previous kernel, at the cost of computing the
    // phase factor once per tile of orders instead of once overall.
    // Consecutive threads take consecutive grid points of the same tile, so
    // they read the same elements of global_vec_w and write contiguous
    // elements of vec_wsum.
    // The simd clause matters on GPUs: GCC maps each OpenMP thread to a
    // whole warp and uses its 32 lanes only for simd loops, so without it
    // 31 of every 32 GPU threads stay idle.
#pragma omp target teams distribute parallel for simd collapse(2)
    for (int jtile = 0; jtile < num_tiles; jtile++) {
      for (long long ixyz = 0; ixyz < nvtot; ixyz++) {
        const int first = jtile * order_tile; // first order of the tile, from jlmf
        const int count = (num_orders - first < order_tile) ? num_orders - first : order_tile;
        const dcomplex *tile_w = global_vec_w + nkvs * first;
        long long iz75 = ixyz / nvxy;
        long long iy70 = (ixyz % nvxy) / nxv;
        long long ix65 = ixyz % nxv;
        double z = _zv[iz75] + frsh;
        double y = _yv[iy70];
        double x = _xv[ix65];
        double sum_re[order_tile], sum_im[order_tile];
        for (int t = 0; t < order_tile; t++) {
          sum_re[t] = 0.0;
          sum_im[t] = 0.0;
        }

        for (long long jy60x55 = 0; jy60x55 < nkvs ; jy60x55++) {
          long long jy60 = jy60x55 / nkv;
          long long jx55 = jy60x55 % nkv;
          long long w_index = (jx55 * nkv) + jy60;
          double vky = vkv[jy60];
          double factor = (jy60 == 0 || jy60 == nkvmo) ? 0.5 : 1.0;
          double vkx, vkz;
          if (jx55 == 0) {
            vkx = -vkv[nkvmo];
            vkz = vec_vkzm[jy60];
            factor *= 0.5;
          } else if (jx55 == nkvmo) {
            vkx = vkv[nkvmo];
            vkz = vec_vkzm[nkvvmo + jy60];
            factor *= 0.5;
          } else {
            vkx = vkv[jx55];
            vkz = vec_vkzm[w_index];
          }
          double angle = vkx * x + vky * y + vkz * z;
          double s, c;
          sincos(angle, &s, &c);
          const double term_re = factor * c;
          const double term_im = factor * s;
          // Complex multiply-add written out in real arithmetic, to avoid
          // the NaN/Inf handling of the C complex product. The loop is fully
          // unrolled so that the sums stay in registers.
#pragma GCC unroll order_tile
          for (int t = 0; t < order_tile; t++) {
            if (t < count) {
              const dcomplex w = tile_w[nkvs * t + w_index];
              const double w_re = __real__ w;
              const double w_im = __imag__ w;
              sum_re[t] += w_re * term_re - w_im * term_im;
              sum_im[t] += w_re * term_im + w_im * term_re;
            }
          }
        } // jy60x55 loop

        for (int t = 0; t < order_tile; t++) {
          if (t < count) {
            long long wsum_index = ((long long)(jlmf - 1 + first + t) * nvtot) + ixyz;
            dcomplex sum;
            __real__ sum = delks_value * sum_re[t];
            __imag__ sum = delks_value * sum_im[t];
            vec_wsum[wsum_index] += sum;
          }
        }
      } // ixyz loop
    } // jtile loop
  } // target region
}
#endif // USE TARGET_OFFLOAD

#ifndef USE_TARGET_OFFLOAD
void cpu_loop(
  dcomplex *vec_wsum, const dcomplex *vec_tt1_wk, double *vkv, double *_xv,
  const int nxv, double *_yv, const int nyv, double *_zv, const int nzv,
  double *vec_vkzm, const int jlmf, const int jlml, const int nkv,
  const int nlmmt, double delks, double frsh
) {
  constexpr int order_tile = FRFME_CPU_ORDER_TILE;
  constexpr int point_block = FRFME_CPU_POINT_BLOCK;
  const long long nvtot = (long long)nxv * nyv * nzv;
  const long long nkvs = (long long)nkv * nkv;
  const long long nkvmo = nkv - 1;
  const long long nvxy = (long long)nxv * nyv;
  const int num_orders = jlml - jlmf + 1;
  const int num_tiles = (num_orders + order_tile - 1) / order_tile;
  const long long num_blocks = (nvtot + point_block - 1) / point_block;

  // Wave vector components and integration weights, in the order of the
  // sum (jy60x55 = jy60 * nkv + jx55, as in offload_loop()).
  double *wave_x = new double[nkvs];
  double *wave_y = new double[nkvs];
  double *wave_z = new double[nkvs];
  double *weight = new double[nkvs];
#pragma omp parallel for
  for (long long jy60x55 = 0; jy60x55 < nkvs; jy60x55++) {
    long long jy60 = jy60x55 / nkv;
    long long jx55 = jy60x55 % nkv;
    double factor = (jy60 == 0 || jy60 == nkvmo) ? 0.5 : 1.0;
    wave_y[jy60x55] = vkv[jy60];
    if (jx55 == 0) {
      wave_x[jy60x55] = -vkv[nkvmo];
      wave_z[jy60x55] = vec_vkzm[jy60];
      factor *= 0.5;
    } else if (jx55 == nkvmo) {
      wave_x[jy60x55] = vkv[nkvmo];
      wave_z[jy60x55] = vec_vkzm[nkvmo * nkv + jy60];
      factor *= 0.5;
    } else {
      wave_x[jy60x55] = vkv[jx55];
      wave_z[jy60x55] = vec_vkzm[jx55 * nkv + jy60];
    }
    weight[jy60x55] = factor;
  }

  // Coefficients rearranged by tile of orders, so that the coefficients of
  // a tile for one wave vector are contiguous. vec_tt1_wk stores the
  // coefficient of order j for the wave vector (jx, jy) at
  // nlmmt * (jy * nkv + jx) + j. Orders past jlml in the last tile are zero.
  const long long tile_size = nkvs * order_tile;
  dcomplex *tile_w = new dcomplex[(long long)num_tiles * tile_size]();
#pragma omp parallel for collapse(2)
  for (int jtile = 0; jtile < num_tiles; jtile++) {
    for (long long jy60x55 = 0; jy60x55 < nkvs; jy60x55++) {
      dcomplex *dest = tile_w + jtile * tile_size + jy60x55 * order_tile;
      for (int t = 0; t < order_tile; t++) {
        int order = jtile * order_tile + t;
        if (order < num_orders) {
          dest[t] = vec_tt1_wk[nlmmt * jy60x55 + jlmf - 1 + order];
        }
      }
    }
  }

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (int jtile = 0; jtile < num_tiles; jtile++) {
    for (long long jblock = 0; jblock < num_blocks; jblock++) {
      const int first_order = jtile * order_tile;
      const int order_count = (num_orders - first_order < order_tile) ? num_orders - first_order : order_tile;
      const long long first_point = jblock * point_block;
      const int point_count = (nvtot - first_point < point_block) ? (int)(nvtot - first_point) : point_block;
      const dcomplex *block_w = tile_w + jtile * tile_size;
      double x[point_block], y[point_block], z[point_block];
      for (int p = 0; p < point_block; p++) {
        // Points past the end of the grid repeat the last point; their sums
        // are not stored.
        long long ixyz = first_point + ((p < point_count) ? p : point_count - 1);
        x[p] = _xv[ixyz % nxv];
        y[p] = _yv[(ixyz % nvxy) / nxv];
        z[p] = _zv[ixyz / nvxy] + frsh;
      }
      double sum_re[order_tile][point_block], sum_im[order_tile][point_block];
      for (int t = 0; t < order_tile; t++) {
        for (int p = 0; p < point_block; p++) {
          sum_re[t][p] = 0.0;
          sum_im[t][p] = 0.0;
        }
      }
      double term_re[point_block], term_im[point_block];
      for (long long jy60x55 = 0; jy60x55 < nkvs; jy60x55++) {
        const double vkx = wave_x[jy60x55];
        const double vky = wave_y[jy60x55];
        const double vkz = wave_z[jy60x55];
        const double factor = weight[jy60x55];
        for (int p = 0; p < point_block; p++) {
          double s, c;
          sincos(vkx * x[p] + vky * y[p] + vkz * z[p], &s, &c);
          term_re[p] = factor * c;
          term_im[p] = factor * s;
        }
        const dcomplex *w = block_w + jy60x55 * order_tile;
        for (int t = 0; t < order_tile; t++) {
          const double w_re = __real__ w[t];
          const double w_im = __imag__ w[t];
#pragma omp simd
          for (int p = 0; p < point_block; p++) {
            sum_re[t][p] += w_re * term_re[p] - w_im * term_im[p];
            sum_im[t][p] += w_re * term_im[p] + w_im * term_re[p];
          }
        }
      } // jy60x55 loop
      for (int t = 0; t < order_count; t++) {
        dcomplex *row = vec_wsum + (long long)(jlmf - 1 + first_order + t) * nvtot + first_point;
        for (int p = 0; p < point_count; p++) {
          dcomplex sum;
          __real__ sum = delks * sum_re[t][p];
          __imag__ sum = delks * sum_im[t][p];
          row[p] = sum;
        }
      }
    } // jblock loop
  } // jtile loop

  delete[] tile_w;
  delete[] wave_x;
  delete[] wave_y;
  delete[] wave_z;
  delete[] weight;
}
#endif // USE_TARGET_OFFLOAD
