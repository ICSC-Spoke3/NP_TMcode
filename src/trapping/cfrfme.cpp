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

#ifdef USE_NVTX
#include <nvtx3/nvToolsExt.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_TARGET_OFFLOAD
#include <cstdlib>

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
 * \param int size_global_vec_w: `int` Size of the global weight vector space.
 * \param vec_tt1_wk: `dcomplex *` Vector of swap object weights.
 * \param int size_vec_tt1_wk: `int` Size of the swap object weights vector.
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

using namespace std;

/*! \brief C++ implementation of FRFME
 *
 *  \param data_file: `string` Name of the input data file.
 *  \param output_path: `string` Directory to write the output files in.
 */
void frfme(string data_file, string output_path) {
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
  string str_target = file_lines[last_read_line++];
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
	message = "ERROR: could not open TEMPTAPE2 file.\n";
	logger.err(message);
      }
    } else {
      message = "ERROR: could not open TFRFME file.\n";
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
    str_target = file_lines[last_read_line++];
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
    str_target = file_lines[last_read_line++];
    for (int cli = 0; cli < 3; cli++) {
      regex_search(str_target, m, re);
      if (cli == 0) an = stod(m.str());
      else if (cli == 1) ff = stod(m.str());
      else if (cli == 2) tra = stod(m.str());
      str_target = m.suffix().str();
    }
    str_target = file_lines[last_read_line++];
    for (int cli = 0; cli < 3; cli++) {
      regex_search(str_target, m, re);
      if (cli == 0) spd = stod(m.str());
      else if (cli == 1) spdfr = stod(m.str());
      else if (cli == 2) exdcl = stod(m.str());
      str_target = m.suffix().str();
    }
    str_target = file_lines[last_read_line++];
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
	throw(UnrecognizedConfigurationException("ERROR: only 'm', 'M', 'e', or 'E' accepted as modes!\n"));
      }
    }
    str_target = m.suffix().str();
    re = regex("[0-9]+");
    regex_search(str_target, m, re);
    int ixi = stoi(m.str());
    string tedf_name = output_path + "/" + namef + ".hd5";
    // Check for run-time options
    bool skip_frfme = false;
    last_read_line++; // skip the end-of-data code
    while (last_read_line < line_count) {
      bool is_parsed = false;
      str_target = file_lines[last_read_line++];
      if (str_target.size() > 12) {
	if (str_target.substr(0, 12).compare("PRECOMPUTED=") == 0) {
	  string precomp_name = output_path + "/" + str_target.substr(12, str_target.size());
	  // TODO: check that the file exists and it is suitable for the calculation.
	  message = "INFO: using precomputed file " + precomp_name + "\n";
	  logger.log(message);
	  skip_frfme = true;
	  is_parsed = true;
	}
      }
      if (!is_parsed) {
	if (str_target.size() > 0) {
	  if (str_target.substr(0, 1).compare("#") != 0) {
	    // ERROR: found an unrecognized option that is not marked as comment
	    throw(UnrecognizedConfigurationException("ERROR: unrecognized option \"" + str_target + "\"!\n"));
	  }
	}
      }
    }
    if (skip_frfme) return;
    ScattererConfiguration *tedf = ScattererConfiguration::from_binary(tedf_name, "HDF5");
#ifdef USE_NVTX
    nvtxRangePush("TEDF data import");
#endif
    if (tedf == NULL) {
      throw(OpenConfigurationFileException("ERROR: could not open " + tedf_name + "hd5!\n"));
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
    message = "INFO: calculating memory requirements...\n";
    logger.log(message);
    swap1_size = Swap1::get_size(lm, nkv);
    size_mb = 1.0 * swap1_size / 1024.0 / 1024.0;
    printf("Swap 1: %.2lg MB\n", size_mb);
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
    size_mb = 1.0 * (swap1_size + swap2_size + tfrfme_size) / 1024.0 / 1024.0;
    sprintf(buffer, "TOTAL: %.2lg MB\n", size_mb);
    message = string(buffer);
    logger.log(message);
    sprintf(buffer, "INFO: nxv = %d, nyv = %d, nzv = %d, nrvc = %d\n", nxv, nyv, nzv, (nxv * nyv * nzv));
    message = string(buffer);
    logger.log(message);
    sprintf(buffer, "INFO: nlmmt = %d, nkv = %d\n", nlmmt, nkv);
    message = string(buffer);
    logger.log(message);
    
    tfrfme = new TFRFME(lmode, lm, nkv, nxv, nyv, nzv);
    if (tfrfme == NULL)
      throw(ObjectAllocationException("ERROR: could not initialize TFRFME structure!\n"));
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
  const int nkvs = nkv * nkv;
  int size_vec_wsum = nlmmt * nrvc;
  int size_global_vec_w = nkvs * (jlml - jlmf + 1);
  int size_vec_tt1_wk = nkvs * nlmmt;
  const dcomplex *vec_tt1_wk = tt1->wk;
  dcomplex *vec_wsum = tfrfme->vec_wsum;
  double *vec_vkzm = tt2->vec_vkzm;
#ifdef USE_TARGET_OFFLOAD
  dcomplex *global_vec_w = (dcomplex *)aligned_alloc(64, size_global_vec_w * sizeof(dcomplex));
#else
  dcomplex *global_vec_w = new dcomplex[size_global_vec_w]();
#endif // USE_TARGET_OFFLOAD
  message = "INFO: looping over " + to_string(jlml - jlmf + 1) + " J iterations.\n";
  logger.log(message);
#ifdef USE_TARGET_OFFLOAD
  t_end = chrono::high_resolution_clock::now();
  elapsed = t_start - t_end;
  frfme_duration = elapsed;
  t_start = chrono::high_resolution_clock::now();
  message = "INFO: computing loop.\n";
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
  sprintf(buffer, "INFO: loop calculation took %lfs.\n", elapsed.count());
  message = string(buffer);
  logger.log(message);
  free(global_vec_w);
#else
  // This code block is compiled if USE_TARGET_OFFLOAD is not defined
#pragma omp parallel for
  for (int j80 = jlmf - 1; j80 < jlml; j80++) {
    dcomplex *vec_w = global_vec_w + nkvs * (j80 - jlmf + 1);
#pragma omp parallel for simd
    for (long long jxy50 = 0; jxy50 < nkvs; jxy50++) {
      long long wk_index = nlmmt * jxy50;
      dcomplex wk_value = vec_tt1_wk[wk_index + j80];
      long long jy50 = jxy50 / nkv;
      long long jx50 = jxy50 % nkv;
      vec_w[(nkv * jx50) + jy50] = wk_value;
    } // jxy50 loop
    long long nvxy = nxv * nyv;
#pragma omp parallel for
    for (long long ixyz = 0; ixyz < nrvc; ixyz++) {
      long long iz75 = ixyz / nvxy;
      long long iy70 = (ixyz % nvxy) / nxv;
      long long ix65 = ixyz % nxv;
      double z = _zv[iz75] + frsh;
      double y = _yv[iy70];
      double x = _xv[ix65];
      dcomplex sumy = cc0;
#pragma omp parallel for simd reduction(+:sumy)
      for (long long jy60x55 = 0; jy60x55 < nkvs ; jy60x55++) {
	long long jy60 = jy60x55 / nkv;
	long long jx55 = jy60x55 % nkv;
	long long w_index = (jx55 * nkv) + jy60;
	double vky = vkv[jy60];
	if (jx55 == 0) {
	  // jx55 = 0: phasf
	  double vkx = vkv[nkv - 1];
	  double vkz = vec_vkzm[jy60];
	  dcomplex phasf = cexp(uim * (-vkx * x + vky * y + vkz * z));
	  dcomplex term = vec_w[jy60] * phasf * 0.5;
	  double factor = (jy60 == 0 || jy60 == nkv - 1) ? 0.5 : 1.0;
	  term *= factor;
	  sumy += term;
	} else if (jx55 == nkv - 1) {
	  // jx55 = nkv - 1: phasl
	  double vkx = vkv[nkv - 1];
	  double vkz = vec_vkzm[(nkv - 1) * nkv + jy60];
	  dcomplex phasl = cexp(uim * (vkx * x + vky * y + vkz * z));
	  dcomplex term = vec_w[(nkv - 1) * nkv + jy60] * phasl * 0.5;
	  double factor = (jy60 == 0 || jy60 == nkv - 1) ? 0.5 : 1.0;
	  term *= factor;
	  sumy += term;
	} else {
	  // 1 <= jx55 < nkv - 1
	  double vkx = vkv[jx55];
	  double vkz = vec_vkzm[(jx55) * nkv + jy60];
	  dcomplex phas = cexp(uim * (vkx * x + vky * y + vkz * z));
	  dcomplex term = vec_w[(jx55) * nkv + jy60] * phas;
	  double factor = (jy60 == 0 || jy60 == nkv - 1) ? 0.5 : 1.0;
	  term *= factor;
	  sumy += term;
	}
      } // jy60x55 loop
      vec_wsum[(j80 * nrvc) + ixyz] = sumy * delks;
    } // ixyz loop
  } // j80 loop
  delete[] global_vec_w;
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
  fprintf(output, " IF JLML < NLMMT, PRESERVE TEMPTAPE1, TEMPTAPE2, AND TFRFRME,\n");
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
  message = "INFO: FRFME took " + to_string(frfme_duration.count()) + "s.\n";
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
  const dcomplex uim = 0.0 + I * 1.0;
  dcomplex cc0 = 0.0 + I * 0.0;

 // Inizializza global_vec_w e vec_wsum sulla CPU in parallelo
#pragma omp parallel for collapse(2)
  for (int j80 = jlmf - 1; j80 < jlml; j80++) {
    for (long long jxy50 = 0; jxy50 < nkvs; jxy50++) {
      int j80_index = j80 - jlmf + 1;
      dcomplex *vec_w = global_vec_w + nkvs * j80_index;
      long long wk_index = nlmmt * jxy50;
      dcomplex wk_value = vec_tt1_wk[wk_index + j80_index];
      long long jy50 = jxy50 / nkv;
      long long jx50 = jxy50 % nkv;
      vec_w[(nkv * jx50) + jy50] = wk_value;
    } // jxy50 loop
  }
  
#pragma omp parallel for
  for (long i = 0; i < size_vec_wsum; i++) {
    vec_wsum[i] = cc0;
  }

#pragma omp target data map(tofrom: vec_wsum[0:size_vec_wsum]) \
  map(to: global_vec_w[0:size_global_vec_w]) \
  map(to: _xv[0:nxv], _yv[0:nyv], _zv[0:nzv]) \
  map(to: vkv[0:nkv], vec_vkzm[0:nkvs])
  {
    // Run the main calculation in a single GPU kernel
    // with parallelization on the largest loops (nvtot & nkvs)
#pragma omp target teams distribute parallel for collapse(2)
    for (long long ixyz = 0; ixyz < nvtot; ixyz++) {
      for (long long jy60x55 = 0; jy60x55 < nkvs ; jy60x55++) {

        long long iz75 = ixyz / nvxy;
        long long iy70 = (ixyz % nvxy) / nxv;
        long long ix65 = ixyz % nxv;
        double z = _zv[iz75] + frsh;
        double y = _yv[iy70];
        double x = _xv[ix65];

        long long jy60 = jy60x55 / nkv;
        long long jx55 = jy60x55 % nkv;
        long long w_index = (jx55 * nkv) + jy60;
        double vky = vkv[jy60];
        double factor = (jy60 == 0 || jy60 == nkvmo) ? 0.5 : 1.0;
        double vkx, vkz;
        dcomplex phas, term;

        if (jx55 == 0) {
          vkx = vkv[nkvmo];
          vkz = vec_vkzm[jy60];
          double angle = -vkx * x + vky * y + vkz * z;
          double s, c;
          sincos(angle, &s, &c);
          phas = c + uim * s;
          term = phas * 0.5;
          term *= factor;
        } else if (jx55 == nkvmo) {
          vkx = vkv[nkvmo];
          vkz = vec_vkzm[nkvvmo + jy60];
          double angle = vkx * x + vky * y + vkz * z;
          double s, c;
          sincos(angle, &s, &c);
          phas = c + uim * s;
          term = phas * 0.5;
          term *= factor;
        } else {
          vkx = vkv[jx55];
          vkz = vec_vkzm[w_index];
          double angle = vkx * x + vky * y + vkz * z;
          double s, c;
          sincos(angle, &s, &c);
          phas = c + uim * s;
          term = phas * factor;
        }
        
	// The last loop can now be serialized, granting for correctness
	// of the results
        for (int j80 = jlmf - 1; j80 < jlml; j80++) {
          int j80_index = j80 - jlmf + 1;
          dcomplex *vec_w = global_vec_w + nkvs * j80_index;
          long long wsum_index = (j80_index * nvtot) + ixyz;
          vec_wsum[wsum_index] += delks * vec_w[w_index] * term;
        }
      } // jy60x55 loop
    } // ixyz loop
  } // target region
}
#endif // USE TARGET_OFFLOAD
