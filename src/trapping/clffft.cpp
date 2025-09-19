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
  string ttms_name = output_path + "/c_TTMS";
  fstream ttms;
  ttms.open(ttms_name, ios::in | ios::binary);
  if (ttms.is_open()) {
    ttms.read(reinterpret_cast<char *>(&is), sizeof(int));
    ttms.read(reinterpret_cast<char *>(&le), sizeof(int));
    ttms.read(reinterpret_cast<char *>(&vks), sizeof(double));
    ttms.read(reinterpret_cast<char *>(&exris), sizeof(double));
    cil->le = le;
    cil->nlem = le * (le + 2);
    cil->nlemt = cil->nlem + cil->nlem;
    if (is >= 2222) { // label 120
      tms = new dcomplex*[le];
      for (int ti = 0; ti < le; ti++) tms[ti] = new dcomplex[3]();
      // QUESTION|WARNING: original code uses LM without defining it. Where does it come from?
      int lm = le;
      for (int i = 0; i < lm; i++) {
	double vreal, vimag;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	tms[i][0] = vreal + vimag * I;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	tms[i][1] = vreal + vimag * I;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	tms[i][2] = vreal + vimag * I;
      } // i loop
    } else if (is >= 1111) { // label 125
      tmsm = new dcomplex[le]();
      tmse = new dcomplex[le]();
      for (int i = 0; i < le; i++) {
	double vreal, vimag;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	tmsm[i] = vreal + vimag * I;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	tmse[i] = vreal + vimag * I;
      } // i loop
    } else if (is >= 0) { // label 135
      am0m = new dcomplex*[cil->nlemt];
      for (int ai = 0; ai < cil->nlemt; ai++) am0m[ai] = new dcomplex[cil->nlemt]();
      for (int i = 0; i < cil->nlemt; i++) {
	for (int j = 0; j < cil->nlemt; j++) {
	  double vreal, vimag;
	  ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	  ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	  am0m[i][j] = vreal + vimag * I;
	} // j loop
      } // i loop
    } else if (is < 0) {
      nvam = le * le + (le * (le + 1) * (le * 2 + 1)) / 3;
      amd = new dcomplex*[nvam];
      for (int ai = 0; ai < nvam; ai++) amd[ai] = new dcomplex[4]();
      for (int i = 0; i < nvam; i++) {
	for (int j = 0; j < 4; j++) {
	  double vreal, vimag;
	  ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	  ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	  amd[i][j] = vreal + vimag * I;
	} // j loop
      } // i loop
      indam = new int*[le];
      int vint;
      for (int ii = 0; ii < le; ii++) indam[ii] = new int[le]();
      for (int i = 0; i < le; i++) {
	for (int j = 0; j < le; j++) {
	  ttms.read(reinterpret_cast<char *>(&vint), sizeof(int));
	  indam[i][j] = vint;
	} // j loop
      } // i loop
      ttms.read(reinterpret_cast<char *>(&vint), sizeof(int));
      cil->mxmpo = vint;
      cil->mxim = vint * 2 - 1;
    }
    // label 150
    ttms.close();
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
	      string tlfff_name = output_path + "/c_TLFFF";
	      tlfff.open(tlfff_name.c_str(), ios::out | ios::binary);
	      if (tlfff.is_open()) {
		tlfff.write(reinterpret_cast<char *>(&lmode), sizeof(int));
		tlfff.write(reinterpret_cast<char *>(&le), sizeof(int));
		tlfff.write(reinterpret_cast<char *>(&nkv), sizeof(int));
		tlfff.write(reinterpret_cast<char *>(&nxv), sizeof(int));
		tlfff.write(reinterpret_cast<char *>(&nyv), sizeof(int));
		tlfff.write(reinterpret_cast<char *>(&nzv), sizeof(int));
		tlfff.write(reinterpret_cast<char *>(&vk), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&exri), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&an), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&ff), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&tra), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&spd), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&frsh), sizeof(double));
		tlfff.write(reinterpret_cast<char *>(&exril), sizeof(double));
		for (int i = 0; i < nxv; i++) {
		  double x = _xv[i];
		  tlfff.write(reinterpret_cast<char *>(&x), sizeof(double));
		}
		for (int i = 0; i < nyv; i++) {
		  double y = _yv[i];
		  tlfff.write(reinterpret_cast<char *>(&y), sizeof(double));
		}
		for (int i = 0; i < nzv; i++) {
		  double z = _zv[i];
		  tlfff.write(reinterpret_cast<char *>(&z), sizeof(double));
		}
		if (jft < 0) goto160 = true;
	      } else { // Should never happen.
		printf("ERROR: could not open TLFFF file.\n");
	      }
	    }
	  }
	  // label 155
	  if (!goto160) {
	    if (jss != 1) {
	      string tlfft_name = output_path + "/c_TLFFT";
	      tlfft.open(tlfft_name.c_str(), ios::out | ios::binary);
	      if (tlfft.is_open()) {
		string outgrid_name = output_path + "/c_grid_scale.txt";
		FILE *outgrid = fopen(outgrid_name.c_str(), "w");
		tlfft.write(reinterpret_cast<char *>(&lmode), sizeof(int));
		tlfft.write(reinterpret_cast<char *>(&le), sizeof(int));
		tlfft.write(reinterpret_cast<char *>(&nkv), sizeof(int));
		tlfft.write(reinterpret_cast<char *>(&nxv), sizeof(int));
		tlfft.write(reinterpret_cast<char *>(&nyv), sizeof(int));
		tlfft.write(reinterpret_cast<char *>(&nzv), sizeof(int));
		tlfft.write(reinterpret_cast<char *>(&vk), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&exri), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&an), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&ff), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&tra), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&spd), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&frsh), sizeof(double));
		tlfft.write(reinterpret_cast<char *>(&exril), sizeof(double));
		for (int i = 0; i < nxv; i++) {
		  double x = _xv[i];
		  if (nxv - i == i + 1) {
		    double ximo = _xv[i - 1];
		    double xipo = _xv[i + 1];
		    if (ximo < 0.0 && xipo > 0.0) {
		      double limo = log10(-ximo);
		      double lipo = log10(xipo);
		      double logi = (x > 0.0) ? log10(x) : log10(-x);
		      if (logi < (limo + lipo) / 2.0 - 10.0) x = 0.0;
		    }
		  }
		  tlfft.write(reinterpret_cast<char *>(&x), sizeof(double));
		  fprintf(outgrid, "  %24.16lE\n", x);
		}
		for (int i = 0; i < nyv; i++) {
		  double y = _yv[i];
		  tlfft.write(reinterpret_cast<char *>(&y), sizeof(double));
		}
		for (int i = 0; i < nzv; i++) {
		  double z = _zv[i];
		  tlfft.write(reinterpret_cast<char *>(&z), sizeof(double));
		}
		fclose(outgrid);
	      } else { // Should never happen.
		printf("ERROR: could not open TLFFT file.\n");
	      }
	    }
	  }
	  // label 160
	  const int nlmm = lm * (lm + 2);
	  const int nlmmt = nlmm + nlmm;
	  const int nrvc = nxv * nyv * nzv;
	  ws = new dcomplex[nlmmt]();
	  if (lm > le) wsl = new dcomplex[nlmmt]();
	  // FORTRAN writes two output formatted files without opening them
	  // explicitly. It is assumed thay can be opened here.
	  string out66_name = output_path + "/c_out66.txt";
	  string out67_name = output_path + "/c_out67.txt";
	  FILE *output66 = fopen(out66_name.c_str(), "w");
	  FILE *output67 = fopen(out67_name.c_str(), "w");
	  for (int iz475 = 0; iz475 < nzv; iz475++) {
	    for (int iy475 = 0; iy475 < nyv; iy475++) {
	      for (int ix475 = 0; ix475 < nxv; ix475++) {
		for (int i = 0; i < nlmmt; i++) {
		  //double vreal, vimag;
		  //binary_input.read(reinterpret_cast<char *>(&vreal), sizeof(double));
		  //binary_input.read(reinterpret_cast<char *>(&vimag), sizeof(double));
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
		    fprintf(
			    output66, " %18.16lE%18.16lE%18.16lE\n",
			    fffe[0], fffs[0], fffe[0] - fffs[0]
			    );
		    fprintf(
			    output66, " %18.16lE%18.16lE%18.16lE\n",
			    fffe[1], fffs[1], fffe[1] - fffs[1]
			    );
		    fprintf(
			    output66, " %18.16lE%18.16lE%18.16lE\n",
			    fffe[2], fffs[2], fffe[2] - fffs[2]
			    );
		  } else { // label 450
		    for (int i = 0; i < 3; i++) {
		      double value = fffe[i] - fffs[i];
		      tlfff.write(reinterpret_cast<char *>(&value), sizeof(double));
		    }
		    if (jtw == 1) {
		      // Writes to 66
		      fprintf(
			      output66, " %5d%4d%4d%15.4lE%15.4lE%15.4lE\n",
			      ix475 + 1, iy475 + 1, iz475 + 1,
			      fffe[0] - fffs[0], fffe[1] - fffs[1], fffe[2] - fffs[2]
			      );
		    }
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
		    fprintf(
			    output67, " %18.16lE%18.16lE%18.16lE\n",
			    ffte[0], ffts[0], ffte[0] - ffts[0]
			    );
		    fprintf(
			    output67, " %18.16lE%18.16lE%18.16lE\n",
			    ffte[1], ffts[1], ffte[1] - ffts[1]
			    );
		    fprintf(
			    output67, " %18.16lE%18.16lE%18.16lE\n",
			    ffte[2], ffts[2], ffte[2] - ffts[2]
			    );
		  } else { // label 470
		    for (int i = 0; i < 3; i++) {
		      double value = ffte[i] - ffts[i];
		      tlfft.write(reinterpret_cast<char *>(&value), sizeof(double));
		    }
		    if (jtw == 1) {
		      // Writes to 67
		      fprintf(
			      output67, " %5d%4d%4d%15.4lE%15.4lE%15.4lE\n",
			      ix475 + 1, iy475 + 1, iz475 + 1,
			      ffte[0] - ffts[0], ffte[1] - ffts[1], ffte[2] - ffts[2]
			      );
		    }
		  }
		  delete[] ffte;
		  delete[] ffts;
		}
	      } // ix475 loop
	    } // iy475 loop
	  } // iz475 loop
	  if (jss != 1) {
	    if (jft <= 0) tlfff.close();
	    if (jft >= 0) tlfft.close();
	  }
	  fclose(output66);
	  fclose(output67);
	}
      }
      delete tfrfme;
    } else {
      printf("ERROR: could not open binary input file %s.\n", binary_name.c_str());
    }
  } else {
    printf("ERROR: could not open TTMS file.\n");
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
    for (int ai = cil->nlemt - 1; ai > -1; ai--) delete[] am0m[ai];
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
