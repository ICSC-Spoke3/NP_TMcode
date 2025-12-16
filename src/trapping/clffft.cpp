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

/*! \file clffft.cpp
 *
 * \brief C++ implementation of LFFFT functions.
 */
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <hdf5.h>
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

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_TRANSITIONMATRIX_H_
#include "../include/TransitionMatrix.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

#ifdef USE_NVTX
#include <nvtx3/nvToolsExt.h>
#endif

using namespace std;

/*! \brief C++ implementation of LFFFT
 *
 *  \param data_file: `string` Name of the input data file.
 *  \param output_path: `string` Directory to write the output files in.
 */
void lffft(string data_file, string output_path) {
#ifdef USE_NVTX
  nvtxRangePush("Running lffft()");
#endif
  chrono::time_point<chrono::high_resolution_clock> t_start = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed;
  const dcomplex uim = 0.0 + 1.0 * I;
  const double sq2i = 1.0 / sqrt(2.0);
  const dcomplex sq2iti = sq2i * uim;

  char buffer[256];
  string message = "INIT";
  Logger logger(LOG_INFO);
  fstream tlfff, tlfft;
  double ****zpv = NULL;
  dcomplex *ac = NULL, *ws = NULL, *wsl = NULL;
  dcomplex *vec_am0m = NULL;
  dcomplex **am0m = NULL;
  dcomplex **amd = NULL;
  int **indam = NULL;
  dcomplex *tmsm = NULL, *tmse = NULL, **tms = NULL;
  int jft, jss, jtw;
  int is, le, nvam = 0;
  double vks, exris;
  CIL *cil = new CIL();
  CCR *ccr = new CCR();

  int num_lines = 0;
  string *file_lines = load_file(data_file, &num_lines);
  regex re = regex("-?[0-9]+");
  smatch m;
  string str_target = file_lines[0];
  for (int mi = 0; mi < 3; mi++) {
    regex_search(str_target, m, re);
    if (mi == 0) jft = stoi(m.str());
    else if (mi == 1) jss = stoi(m.str());
    else if (mi == 2) jtw = stoi(m.str());
    str_target = m.suffix().str();
  } // mi loop
  int last_read_line = 2; // skip the end-of-data code
  double host_ram_gb = 0.0;
  // Parse runtime options.
  while (last_read_line < num_lines) {
    bool is_parsed = false;
    string str_target = file_lines[last_read_line++];
    if (str_target.size() > 12) {
      if (str_target.substr(0, 12).compare("HOST_RAM_GB=") == 0) {
	host_ram_gb = (double)stod(str_target.substr(12, str_target.length()));
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
  string ttms_name = output_path + "/c_TTMS.hd5";
  // fstream ttms;
  // ttms.open(ttms_name, ios::in | ios::binary);
  // if (ttms.is_open()) {
  TransitionMatrix *ttms = TransitionMatrix::from_binary(ttms_name, "HDF5");
  if (ttms == NULL)
    throw(ObjectAllocationException("ERROR: could not load the T-matrix file!\n"));
  // ttms.read(reinterpret_cast<char *>(&is), sizeof(int));
  // ttms.read(reinterpret_cast<char *>(&le), sizeof(int));
  // ttms.read(reinterpret_cast<char *>(&vks), sizeof(double));
  // ttms.read(reinterpret_cast<char *>(&exris), sizeof(double));
  is = ttms->is;
  le = ttms->l_max;
  vks = ttms->vk;
  exris = ttms->exri;
  cil->le = le;
  cil->nlem = le * (le + 2);
  cil->nlemt = cil->nlem + cil->nlem;
  if (is >= 2222) { // label 120
    // READING OF T-MATRIX WITH IS = 2222 IS NOT SUPPORTED!
    // tms = new dcomplex*[le];
    // for (int ti = 0; ti < le; ti++) tms[ti] = new dcomplex[3]();
    // int lm = le;
    // for (int i = 0; i < lm; i++) {
    //   double vreal, vimag;
    //   ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
    //   ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
    //   tms[i][0] = vreal + vimag * I;
    //   ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
    //   ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
    //   tms[i][1] = vreal + vimag * I;
    //   ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
    //   ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
    //   tms[i][2] = vreal + vimag * I;
    //} // i loop
    throw(UnrecognizedConfigurationException("ERROR: T-matrix with IS>=2222 not supported!\n"));
  } else if (is >= 1111) { // label 125
    tmsm = new dcomplex[le]();
    tmse = new dcomplex[le]();
    for (int i = 0; i < le; i++) {
      tmsm[i] = ttms->elements[2 * i];
      tmse[i] = ttms->elements[2 * i + 1];
    } // i loop
  } else if (is >= 0) { // label 135
    vec_am0m = new dcomplex[cil->nlemt * cil->nlemt];
    am0m = new dcomplex*[cil->nlemt];
    for (int ai = 0; ai < cil->nlemt; ai++) am0m[ai] = vec_am0m + cil->nlemt * ai;
    for (int i = 0; i < cil->nlemt; i++) {
      for (int j = 0; j < cil->nlemt; j++) {
	am0m[i][j] = ttms->elements[cil->nlemt * i + j];
      } // j loop
    } // i loop
  } else if (is < 0) {
    // nvam = le * le + (le * (le + 1) * (le * 2 + 1)) / 3;
    // amd = new dcomplex*[nvam];
    // for (int ai = 0; ai < nvam; ai++) amd[ai] = new dcomplex[4]();
    // for (int i = 0; i < nvam; i++) {
    //   for (int j = 0; j < 4; j++) {
    // 	double vreal, vimag;
    // 	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
    // 	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
    // 	amd[i][j] = vreal + vimag * I;
    //   } // j loop
    // } // i loop
    // indam = new int*[le];
    // int vint;
    // for (int ii = 0; ii < le; ii++) indam[ii] = new int[le]();
    // for (int i = 0; i < le; i++) {
    //   for (int j = 0; j < le; j++) {
    // 	ttms.read(reinterpret_cast<char *>(&vint), sizeof(int));
    // 	indam[i][j] = vint;
    //   } // j loop
    // } // i loop
    // ttms.read(reinterpret_cast<char *>(&vint), sizeof(int));
    // cil->mxmpo = vint;
    // cil->mxim = vint * 2 - 1;
    throw(UnrecognizedConfigurationException("ERROR: T-matrix with IS<0 not supported!\n"));
  }
  // label 150
  // ttms.close();
  delete ttms;
  TFRFME *tfrfme = NULL;
  string binary_name;
  if (jss != 1) binary_name = output_path + "/c_TFRFME.hd5";
  else binary_name = output_path + "/c_TWS";
  tfrfme = TFRFME::from_binary(binary_name, "HDF5");
  if (tfrfme != NULL) {
    int lmode, lm, nkv, nxv, nyv, nzv;
    double vk, exri, an, ff, tra;
    double spd, frsh, exril;
    lmode = tfrfme->lmode;
    lm = tfrfme->lm;
    nkv = tfrfme->nkv;
    nxv = tfrfme->nxv;
    nyv = tfrfme->nyv;
    nzv = tfrfme->nzv;
    double *_xv = tfrfme->get_x();
    double *_yv = tfrfme->get_y();
    double *_zv = tfrfme->get_z();
    if (host_ram_gb > 0.0) {
      double output_size_gb = TrappingOutputInfo::get_size(jft, jss, nxv, nyv, nzv)
	/ 1024.0 / 1024.0 / 1024.0;
      if (output_size_gb >= host_ram_gb) {
	throw(ObjectAllocationException("Output data do not fit in host RAM!\n"));
      }
    }
    TrappingOutputInfo *toi = new TrappingOutputInfo(jft, jss, jtw, nxv, nyv, nzv);
    if (lm >= le) {
      vk = tfrfme->vk;
      exri = tfrfme->exri;
      an = tfrfme->an;
      ff = tfrfme->ff;
      tra = tfrfme->tra;
      if (vk == vks && exri == exris) {
	spd = tfrfme->spd;
	frsh = tfrfme->frsh;
	exril = tfrfme->exril;
	bool goto160 = false;
	if (jft <= 0) {
	  zpv = new double***[le];
	  for (int zi = 0; zi < le; zi++) {
	    zpv[zi] = new double**[3];
	    for (int zj = 0; zj < 3; zj++) {
	      zpv[zi][zj] = new double*[2];
	      for (int zk = 0; zk < 2; zk++) zpv[zi][zj][zk] = new double[2]();
	    } // zj loop
	  } // zi loop
	  thdps(le, zpv);
	  double exdc = exri * exri;
	  double sqk = vk * vk * exdc;
	  ccr->cof = 1.0 / sqk;
	  ccr->cimu = ccr->cof / sqrt(2.0);
	  if (jss != 1) {
	    toi->set_param("lmode", 1.0 * lmode);
	    toi->set_param("le", 1.0 * le);
	    toi->set_param("nkv", 1.0 * nkv);
	    toi->set_param("vk", vk);
	    toi->set_param("exri", exri);
	    toi->set_param("an", an);
	    toi->set_param("ff", ff);
	    toi->set_param("tra", tra);
	    toi->set_param("spd", spd);
	    toi->set_param("frsh", frsh);
	    toi->set_param("exril", exril);
	    for (int i = 0; i < nxv; i++) {
	      toi->vec_x[i] = _xv[i];
	    }
	    for (int i = 0; i < nyv; i++) {
	      toi->vec_y[i] = _yv[i];
	    }
	    for (int i = 0; i < nzv; i++) {
	      toi->vec_z[i] = _zv[i];
	    }
	      if (jft < 0) goto160 = true;
	  }
	}
	// label 155
	// label 160
	const int nlmm = lm * (lm + 2);
	const int nlmmt = nlmm + nlmm;
	const int nrvc = nxv * nyv * nzv;
	ws = new dcomplex[nlmmt]();
	if (lm > le) wsl = new dcomplex[nlmmt]();
	for (int iz475 = 0; iz475 < nzv; iz475++) {
	  for (int iy475 = 0; iy475 < nyv; iy475++) {
	    for (int ix475 = 0; ix475 < nxv; ix475++) {
	      int out_index = nyv * nxv * iz475 + nxv * iy475 + ix475;
	      for (int i = 0; i < nlmmt; i++) {
		int row = i;
		int col = (nyv * nxv * iz475) + (nxv * iy475) + ix475;
		dcomplex value = tfrfme->vec_wsum[nrvc * row + col];
		if (lm <= le) {
		  ws[i] = value;
		} else { // label 170
		  wsl[i] = value;
		  for (int i175 = 0; i175 < cil->nlem; i175++) {
		    int ie = i175 + cil->nlem;
		    int iel = i175 + nlmm;
		    ws[i175] = wsl[i175];
		    ws[ie] = wsl[iel];
		  } // i175 loop
		}
	      } // i loop
	      // label 180
	      if (is != 2222) {
		if (is != 1111) {
		  if (is > 0) { // Goes to 305
		    if (ac) delete[] ac;
		    ac = new dcomplex[cil->nlemt]();
		    camp(ac, am0m, ws, cil);
		    // Goes to 445
		  } else if (is < 0) { // Goes to 405
		    if (ac) delete[] ac;
		    ac = new dcomplex[cil->nlemt]();
		    czamp(ac, amd, indam, ws, cil);
		    // Goes to 445
		  }
		} else {
		  if (ac) delete[] ac;
		  ac = new dcomplex[cil->nlemt]();
		  samp(ac, tmsm, tmse, ws, cil);
		  // Goes to 445
		}
	      } else {
		ac = new dcomplex[cil->nlemt]();
		sampoa(ac, tms, ws, cil);
		// Goes to 445
	      }
	      bool goto475 = false;
	      // label 445
	      if (jft <= 0) {
		double *fffe = new double[3]();
		double *fffs = new double[3]();
		ffrf(zpv, ac, ws, fffe, fffs, cil, ccr);
		if (jss == 1) {
		  // Writes to 66
		  toi->vec_csf1[3 * out_index] = fffe[0];
		  toi->vec_csf2[3 * out_index] = fffs[0];
		  toi->vec_csf3[3 * out_index] = fffe[0] - fffs[0];
		  toi->vec_csf1[3 * out_index + 1] = fffe[1];
		  toi->vec_csf2[3 * out_index + 1] = fffs[1];
		  toi->vec_csf3[3 * out_index + 1] = fffe[1] - fffs[1];
		  toi->vec_csf1[3 * out_index + 2] = fffe[2];
		  toi->vec_csf2[3 * out_index + 2] = fffs[2];
		  toi->vec_csf3[3 * out_index + 2] = fffe[2] - fffs[2];
		} else { // label 450, jss != 1
		  toi->vec_csf1[out_index] = fffe[0] - fffs[0];
		  toi->vec_csf2[out_index] = fffe[1] - fffs[1];
		  toi->vec_csf3[out_index] = fffe[2] - fffs[2];
		}
		if (jft < 0) goto475 = true;
		delete[] fffe;
		delete[] fffs;
	      }
	      // label 460
	      if (!goto475) {
		double *ffte = new double[3]();
		double *ffts = new double[3]();
		ffrt(ac, ws, ffte, ffts, cil);
		if (jss == 1) {
		  // Writes to 67
		  toi->vec_cst1[3 * out_index] = ffte[0];
		  toi->vec_cst2[3 * out_index] = ffts[0];
		  toi->vec_cst3[3 * out_index] = ffte[0] - ffts[0];
		  toi->vec_cst1[3 * out_index + 1] = ffte[1];
		  toi->vec_cst2[3 * out_index + 1] = ffts[1];
		  toi->vec_cst3[3 * out_index + 1] = ffte[1] - ffts[1];
		  toi->vec_cst1[3 * out_index + 2] = ffte[2];
		  toi->vec_cst2[3 * out_index + 2] = ffts[2];
		  toi->vec_cst3[3 * out_index + 2] = ffte[2] - ffts[2];
		} else { // label 470, jss != 1
		  toi->vec_cst1[out_index] = ffte[0] - ffts[0];
		  toi->vec_cst2[out_index] = ffte[1] - ffts[1];
		  toi->vec_cst3[out_index] = ffte[2] - ffts[2];
		}
		delete[] ffte;
		delete[] ffts;
	      }
	    } // ix475 loop
	  } // iy475 loop
	} // iz475 loop
      } else { // vk != vks || exri != exris
	message = "ERROR: T-matrix mismatch!\n";
	if (vk != vks) message = "ERROR: calculation wavelength not matching T-matrix wavelength!\n";
	else if (exri != exris) message = "ERROR: external medium not matching T-matrix external medium!\n";
	throw(UnrecognizedConfigurationException(message));
      }
    } else { // lm < le: ERROR!
      message = "ERROR: calculation order " + to_string(lm) + " smaller than T-matrix order "
	+ to_string(le) + "!\n";
      throw(UnrecognizedConfigurationException(message));
    }
    toi->write(output_path + "/c_LFFFT.hd5", "HDF5");
    if (jtw == 1) toi->write(output_path + "/c_", "ASCII");
    delete tfrfme;
    delete toi;
  } else {
    printf("ERROR: could not open binary input file %s.\n", binary_name.c_str());
  }
  
  // Clean up memory
  if (ac != NULL) delete[] ac;
  if (ws != NULL) delete[] ws;
  if (wsl != NULL) delete[] wsl;
  if (tmsm != NULL) delete[] tmsm;
  if (tmse != NULL) delete[] tmse;
  if (tms != NULL) {
    for (int ti = le - 1; ti > -1; ti--) delete[] tms[ti];
    delete[] tms;
  }
  if (am0m != NULL) {
    //for (int ai = cil->nlemt - 1; ai > -1; ai--) delete[] am0m[ai];
    delete[] vec_am0m;
    delete[] am0m;
  }
  if (amd != NULL) {
    for (int ai = nvam - 1; ai > -1; ai--) delete[] amd[ai];
    delete[] amd;
  }
  if (indam != NULL) {
    for (int ii = le - 1; ii > -1; ii--) delete[] indam[ii];
    delete[] indam;
  }
  if (zpv != NULL) {
    for (int zi = le - 1; zi > -1; zi--) {
      for (int zj = 2; zj > -1; zj--) {
	for (int zk = 1; zk > -1; zk--) delete[] zpv[zi][zj][zk];
	delete[] zpv[zi][zj];
      } // zj loop
      delete[] zpv[zi];
    } // zi loop
    delete[] zpv;
  }
  delete cil;
  delete ccr;
  delete[] file_lines;
  elapsed = chrono::high_resolution_clock::now() - t_start;
  message = "INFO: LFFT took " + to_string(elapsed.count()) + "s.\n";
  logger.log(message);
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}
