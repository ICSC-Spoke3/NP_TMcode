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

/*! \file tra_subs.cpp
 *
 * \brief C++ implementation of TRAPPING subroutines.
 */
#include <cmath>
#include <fstream>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
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

#ifdef _OPENMP
#include <omp.h>
#endif

void camp(dcomplex *ac, dcomplex **am0m, dcomplex *ws, CIL *cil) {
  for (int j = 0; j < cil->nlemt; j++) {
    for (int i = 0; i < cil->nlemt; i++) {
      ac[j] += (am0m[j][i] * ws[i]);
    } // i loop
  } // j loop
}

void czamp(dcomplex *ac, dcomplex **amd, int **indam, dcomplex *ws, CIL *cil) {
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex summ, sume;
  for (int im20 = 1; im20 <= cil->mxim; im20++) {
    int m = im20 - cil->mxmpo;
    int abs_m = (m < 0) ? -m : m;
    int lmn = (abs_m > 1) ? abs_m : 1;
    for (int l20 = lmn; l20 <= cil->le; l20++) {
      int i = m + l20 * (l20 + 1);
      int ie = i + cil->nlem;
      summ = cc0;
      sume = cc0;
      for (int ls15 = lmn; ls15 <= cil->le; ls15++) {
	int is = m + ls15 * (ls15 + 15) - 1;
	int ise = is + cil->nlem;
	int num = indam[l20 - 1][ls15 - 1] + m - 1;
	summ += (amd[num][0] * ws[is] + amd[num][1] * ws[ise]);
	sume += (amd[num][2] * ws[is] + amd[num][3] * ws[ise]);
      } // ls15 loop
      ac[i - 1] = summ;
      ac[ie - 1] = sume;
    } // l20 loop
  } // im20 loop
}

void ffrf(
	  double ****zpv, dcomplex *ac, dcomplex *ws, double *fffe,
	  double *fffs, CIL *cil, CCR *ccr
) {
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex uimmp, summ, sume, suem, suee;
  dcomplex *gap = new dcomplex[3]();

  for (int imu50 = 1; imu50 <= 3; imu50++) {
    int mu = imu50 - 2;
    gap[imu50 - 1] = cc0;
    for (int l40 = 1; l40 <= cil->le; l40++) {
      int lpo = l40 + 1;
      int ltpo = lpo + l40;
      int imm = l40 * lpo;
      for (int ilmp40 = 1; ilmp40 <= 3; ilmp40++) {
	if ((l40 == 1 && ilmp40 == 1) || (l40 == cil->le && ilmp40 == 3)) continue; // ilmp40 loop
	int lmpml = ilmp40 - 2;
	int lmp = l40 + lmpml;
	uimmp = uim * (-1.0 * lmpml);
	int impmmmp = lmp * (lmp + 1);
	for (int im30 = 1; im30 <= ltpo; im30++) {
	  int m = im30 - lpo;
	  int mmp = m - mu;
	  int abs_mmp = (mmp < 0) ? -mmp : mmp;
	  if (abs_mmp <= lmp) {
	    int i = imm + m;
	    int ie = i + cil->nlem;
	    int imp = impmmmp + mmp;
	    int impe = imp + cil->nlem;
	    double cgc = cg1(lmpml, mu, l40, m);
	    summ = dconjg(ws[i - 1]) * ac[imp - 1];
	    sume = dconjg(ws[i - 1]) * ac[impe - 1];
	    suem = dconjg(ws[ie - 1]) * ac[imp - 1];
	    suee = dconjg(ws[ie - 1]) * ac[impe - 1];
	    if (lmpml != 0) {
	      summ *= uimmp;
	      sume *= uimmp;
	      suem *= uimmp;
	      suee *= uimmp;
	    }
	    // label 25
	    gap[imu50 - 1] += (cgc * (
				      summ * zpv[l40 - 1][ilmp40 - 1][0][0] +
				      sume * zpv[l40 - 1][ilmp40 - 1][0][1] +
				      suem * zpv[l40 - 1][ilmp40 - 1][1][0] +
				      suee * zpv[l40 - 1][ilmp40 - 1][1][1]
				      )
			       );
	  }
	} // im30 loop
      } // ilmp40
    } // l40 loop
  } // imu50 loop
  sume = -gap[0] * ccr->cimu;
  suee = -gap[1] * ccr->cof;
  suem = -gap[2] * ccr->cimu;
  fffe[0] = real(sume - suem);
  fffe[1] = real((sume + suem) * uim);
  fffe[2] = real(suee);

  for (int imu90 = 1; imu90 <= 3; imu90++) {
    int mu = imu90 - 2;
    gap[imu90 - 1] = cc0;
    for (int l80 = 1; l80 <= cil->le; l80++) {
      int lpo = l80 + 1;
      int ltpo = lpo + l80;
      int imm = l80 * lpo;
      for (int ilmp80 = 1; ilmp80 <= 3; ilmp80++) {
	if ((l80 == 1 && ilmp80 == 1) || (l80 == cil->le && ilmp80 == 3)) continue; // ilmp80 loop
	int lmpml = ilmp80 - 2;
	int lmp = l80 + lmpml;
	uimmp = uim * (-1.0 * lmpml);
	int impmmmp = lmp * (lmp + 1);
	for (int im70 = 1; im70 <= ltpo; im70++) {
	  int m = im70 - lpo;
	  int mmp = m - mu;
	  int abs_mmp = (mmp < 0) ? -mmp : mmp;
	  if (abs_mmp <= lmp) {
	    int i = imm + m;
	    int ie = i + cil->nlem;
	    int imp = impmmmp + mmp;
	    int impe = imp + cil->nlem;
	    double cgc = cg1(lmpml, mu, l80, m);
	    summ = dconjg(ac[i - 1]) * ac[imp - 1];
	    sume = dconjg(ac[i - 1]) * ac[impe - 1];
	    suem = dconjg(ac[ie - 1]) * ac[imp - 1];
	    suee = dconjg(ac[ie - 1]) * ac[impe - 1];
	    if (lmpml != 0) {
	      summ *= uimmp;
	      sume *= uimmp;
	      suem *= uimmp;
	      suee *= uimmp;
	    }
	    // label 65
	    gap[imu90 - 1] += (cgc * (
				      summ * zpv[l80 - 1][ilmp80 - 1][0][0] +
				      sume * zpv[l80 - 1][ilmp80 - 1][0][1] +
				      suem * zpv[l80 - 1][ilmp80 - 1][1][0] +
				      suee * zpv[l80 - 1][ilmp80 - 1][1][1]
				      )
			       );
	  }
	} // im70 loop
      } // ilmp80 loop
    } // l80 loop
  } // imu90 loop
  sume = gap[0] * ccr->cimu;
  suee = gap[1] * ccr->cof;
  suem = gap[2] * ccr->cimu;
  fffs[0] = real(sume - suem);
  fffs[1] = real((sume + suem) * uim);
  fffs[2] = real(suee);
  delete[] gap;
}

void ffrt(dcomplex *ac, dcomplex *ws, double *ffte, double *ffts, CIL *cil) {
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  const double sq2i = 1.0 / sqrt(2.0);
  const dcomplex sq2iti = uim * sq2i;
  dcomplex aca, acw;
  dcomplex *ctqce, *ctqcs;

  ctqce = new dcomplex[3]();
  ctqcs = new dcomplex[3]();
  for (int l60 = 1; l60 <= cil->le; l60++) {
    int lpo = l60 + 1;
    int il = l60 * lpo;
    int ltpo = l60 + lpo;
    for (int im60 = 1; im60 <= ltpo; im60++) {
      double rmu;
      int m = im60 - lpo;
      int i = m + il;
      int ie = i + cil->nlem;
      int mmmu = m + 1;
      int mmmmu = (mmmu < 0) ? -mmmu: mmmu;
      if (mmmmu <= l60) {
	int immu = mmmu + il;
	int immue = immu + cil->nlem;
	rmu = -sqrt(1.0 * (l60 + mmmu) * (l60 - m)) * sq2i;
	acw = dconjg(ac[i - 1]) * ws[immu - 1] + dconjg(ac[ie - 1]) * ws[immue - 1];
	aca = dconjg(ac[i - 1]) * ac[immu - 1] + dconjg(ac[ie - 1]) * ac[immue - 1];
	ctqce[0] += (acw * rmu);
	ctqcs[0] += (aca * rmu);
      }
      // label 30
      rmu = -1.0 * m;
      acw = dconjg(ac[i - 1]) * ws[i - 1] + dconjg(ac[ie - 1]) * ws[ie - 1];
      aca = dconjg(ac[i - 1]) * ac[i - 1] + dconjg(ac[ie - 1]) * ac[ie - 1];
      ctqce[1] += (acw * rmu);
      ctqcs[1] += (aca * rmu);
      mmmu = m - 1;
      mmmmu = (mmmu < 0) ? -mmmu: mmmu;
      if (mmmmu <= l60) {
	int immu = mmmu + il;
	int immue = immu + cil->nlem;
	rmu = sqrt(1.0 * (l60 - mmmu) * (l60 + m)) * sq2i;
	acw = dconjg(ac[i - 1]) * ws[immu - 1] + dconjg(ac[ie - 1]) * ws[immue - 1];
	aca = dconjg(ac[i - 1]) * ac[immu - 1] + dconjg(ac[ie - 1]) * ac[immue - 1];
	ctqce[2] += (acw * rmu);
	ctqcs[2] += (aca * rmu);
      }
    } // im60 loop
  } // l60 loop
  ffte[0] = real(ctqce[0] - ctqce[2]) * sq2i;
  ffte[1] = real(sq2iti * (ctqce[0] + ctqce[2]));
  ffte[2] = real(ctqce[1]);
  ffts[0] = -sq2i * real(ctqcs[0] - ctqcs[2]);
  ffts[1] = -1.0 * real(sq2iti * (ctqcs[0] + ctqcs[2]));
  ffts[2] = -1.0 * real(ctqcs[1]);

  delete[] ctqce;
  delete[] ctqcs;
}

dcomplex *frfmer(
  int nkv, double vkm, double vknmx, double apfafa, double tra,
  double spd, double rir, double ftcn, int le, int lmode, double pmf,
  Swap1 *tt1, Swap2 *tt2
) {
  const int nlemt = le * (le + 2) * 2;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  double sq = vkm * vkm;
  double *_vkv = tt2->get_vector();
  double *vec_vkzm = tt2->vec_vkzm;
  dcomplex *wk = new dcomplex[nlemt]();
  for (int jy90 = 0; jy90 < nkv; jy90++) {
    double vky = _vkv[jy90];
    double sqy = vky * vky;
    for (int jx80 = 0; jx80 < nkv; jx80++) {
      double vkx = _vkv[jx80];
      double sqx = vkx * vkx;
      double sqn = sqx + sqy;
      double vkn = sqrt(sqn);
      if (vkn <= vknmx) {
	double vkz = sqrt(sq - sqn);
	wamff(wk, vkx, vky, vkz, le, apfafa, tra, spd, rir, ftcn, lmode, pmf);
	for (int j = 0; j < nlemt; j++) tt1->append(wk[j]);
	vec_vkzm[nkv * jx80 + jy90] = vkz;
      } else { // label 50
	for (int j = 0; j < nlemt; j++) tt1->append(cc0);
	vec_vkzm[nkv * jx80 + jy90] = 0.0;
      }
    } // jx80 loop
  } // jy90 loop
  return wk;
}

void pwmalp(dcomplex **w, double *up, double *un, dcomplex *ylm, int lw) {
  dcomplex cp1, cm1, cp2, cm2, cl;
  const dcomplex uim = 0.0 + 1.0 * I;
  const double four_pi = acos(0.0) * 8.0;
  const int nlwm = lw * (lw + 2);
  cm1 = 0.5 * (up[0] + up[1] * I);
  cp1 = 0.5 * (up[0] - up[1] * I);
  double cz1 = up[2];
  cm2 = 0.5 * (un[0] + un[1] * I);
  cp2 = 0.5 * (un[0] - un[1] * I);
  double cz2 =un[2];
  for (int l20 = 1; l20 <= lw; l20++) {
    int lf = l20 + 1;
    int lftl = lf * l20;
    double x = 1.0 * lftl;
    dcomplex cl = (four_pi / sqrt(x)) * cpow(uim, 1.0 * l20);
    int mv = l20 + lf;
    int m = -lf;
    for (int mf20 = 1; mf20 <= mv; mf20++) {
      m++;
      int k = lftl + m;
      x = 1.0 * (lftl - m * (m + 1));
      double cp = sqrt(x);
      x = 1.0 * (lftl - m * (m - 1));
      double cm = sqrt(x);
      double cz = 1.0 * m;
      w[k - 1][0] = dconjg(cp1 * cp * ylm[k + 1] + cm1 * cm * ylm[k - 1] + cz1 * cz * ylm[k]) * cl;
      w[k - 1][1] = dconjg(cp2 * cp * ylm[k + 1] + cm2 * cm * ylm[k - 1] + cz2 * cz * ylm[k]) * cl;
    } // mf20 loop
  } // l20 loop
  for (int k30 = 0; k30 < nlwm; k30++) {
    int i = k30 + nlwm;
    w[i][0] = uim * w[k30][1];
    w[i][1] = -uim * w[k30][0];
  } // k30 loop
}

void samp(dcomplex *ac, dcomplex *tmsm, dcomplex *tmse, dcomplex *ws, CIL *cil) {
  int i = 0;
  for (int l20 = 0; l20 < cil->le; l20++) {
    int l = l20 + 1;
    int ltpo = l + l + 1;
    for (int im20 = 0; im20 < ltpo; im20++) {
      int ie = i + cil->nlem;
      ac[i] = tmsm[l20] * ws[i];
      ac[ie] = tmse[l20] * ws[ie];
      i++;
    } // im20 loop
  } // l20 loop
}

void sampoa(dcomplex *ac, dcomplex **tms, dcomplex *ws, CIL *cil) {
  dcomplex **tm = new dcomplex*[2];
  tm[0] = new dcomplex[2]();
  tm[1] = new dcomplex[2]();
  int i = 0;
  for (int l20 = 0; l20 < cil->le; l20++) {
    tm[0][0] = tms[l20][0];
    tm[0][1] = tms[l20][1];
    tm[1][1] = tms[l20][2];
    tm[1][0] = tm[0][1];
    int l = l20 + 1;
    int ltpo = l + l + 1;
    for (int im20 = 0; im20 < ltpo; im20++) {
      int ie = i + cil->nlem;
      ac[i] = tm[0][0] * ws[i] + tm[0][1] * ws[ie];
      ac[ie] = tm[1][0] * ws[i] + tm[1][1] * ws[ie];
      i++;
    } // im20 loop
  } // l20 loop
  delete[] tm[1];
  delete[] tm[0];
  delete[] tm;
}

void wamff(
	   dcomplex *wk, double x, double y, double z, int lm, double apfafa,
	   double tra, double spd, double rir, double ftcn, int lmode, double pmf
) {
  const int nlmm = lm * (lm + 2);
  const int nlmmt = 2 * nlmm;
  const int nlmp = nlmm + 2;
  dcomplex **w, *ylm;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex cfam, cf1, cf2;
  double rho, cph, sph, cth, sth, r;
  double ftc1, ftc2;
  double *up = new double[3];
  double *un = new double[3];
  w = new dcomplex*[nlmmt];
  for (int wi = 0; wi < nlmmt; wi++) w[wi] = new dcomplex[2]();
  ylm = new dcomplex[nlmp]();
  bool onx = (y == 0.0);
  bool ony = (x == 0.0);
  bool onz = (onx && ony);
  if (!onz) {
    if (!onx) {
      if (!ony) {
	rho = sqrt(x * x + y * y);
	cph = x / rho;
	sph = y / rho;
	// goes to 15
      } else { // label 13
	rho = (y < 0.0) ? -y : y;
	cph = 0.0;
	sph = (y < 0.0) ? -1.0 : 1.0;
	// goes to 15
      }
    } else { // label 12
      rho = (x < 0.0) ? -x : x;
      cph = (x < 0.0) ? -1.0 : 1.0;
      sph = 0.0;
      // goes to 15
    }
  } else { // label 10
    cph = 1.0;
    sph = 0.0;
    // goes to 15
  }
  // label 15
  if (z == 0.0) {
    if (!onz) {
      r = rho;
      cth = 0.0;
      sth = 1.0;
    } else { // label 17
      r = 0.0;
      cth = 1.0;
      sth = 0.0;
    }
  } else { // label 18
    if (!onz) {
      r = sqrt(rho * rho + z * z);
      cth = z / r;
      sth = rho / r;
    } else { // label 20
      r = sqrt(z * z);
      cth = (z < 0.0)? -1.0: 1.0;
      sth = 0.0;
    }
  }
  if (lmode == 0 || sth != 0.0) { // label 25
    bool skip62 = false;
    ylm[nlmp - 1] = cc0;
    sphar(cth, sth, cph, sph, lm, ylm);
    up[0] = cph * cth;
    up[1] = sph * cth;
    up[2] = -sth;
    un[0] = -sph;
    un[1] = cph;
    un[2] = 0.0;
    pwmalp(w, up, un, ylm, lm);
    double apfa = sth * apfafa;
    if (spd > 0.0) {
      double sthl = sth * rir;
      double cthl = sqrt(1.0 - sthl * sthl);
      double arg = spd * (z - (r / rir) * cthl);
      cfam = (tra * cexp(-apfa * apfa) / sqrt(cthl)) * cexp(uim * arg);
      if (lmode == 0) {
	if (sth == 0.0) { // label 45
	  ftc1 = ftcn;
	  ftc2 = ftcn;
	  // goes to 48
	}
      } else if (lmode == 1) { // label 46
	cfam *= ((cph + uim * sph) * sth * pmf);
	ftc1 = 2.0 * cthl / (cthl * rir + cth);
	ftc2 = 2.0 * cthl / (cthl + cth * rir);
	// follows on to 48
      } else if (lmode == 2) { // label 50
	ftc1 = 2.0 * cthl / (cthl * rir + cth);
	cfam *= (sth * pmf * ftc1);
	for (int i52 = 0; i52 < nlmmt; i52++) wk[i52] = w[i52][0] * cfam;
	// returns
	skip62 = true;
      } else if (lmode == 3) { // label 53
	ftc2 = 2.0 * cthl / (cthl  + cth * rir);
	cfam *= (sth * pmf * ftc2);
	for (int i55 = 0; i55 < nlmmt; i55++) wk[i55] = w[i55][1] * cfam;
	// returns
	skip62 = true;
      }
      if (lmode == 0 || lmode == 1) { //label 48
	cf1 = cph * ftc1 * cfam;
	cf2 = -sph * ftc2 * cfam;
	// goes to 62
	skip62 = false;
      }
    } else { // label 57
      double fam = tra * exp(-apfa * apfa) / sqrt(cth);
      if (lmode == 0) {
	double f1 = cph * fam;
	double f2 = -sph * fam;
	for (int i58 = 0; i58 < nlmmt; i58++) wk[i58] = w[i58][0] * f1 + w[i58][1] * f2;
	// returns
	skip62 = true;
      } else if (lmode == 1) { // label 60
	cfam = (pmf * sth * fam) * (cph * uim * sph);
	cf1 = cph * cfam;
	cf2 = -sph * cfam;
	// follows on to 62
	skip62 = false;
      } else if (lmode == 2) { // label 65
	fam *= (pmf * sth);
	for (int i67 = 0; i67 < nlmmt; i67++) wk[i67] = w[i67][0] * fam;
	// returns
	skip62 = true;
      } else if (lmode == 3) { // label 68
	fam *= (pmf * sth);
	for (int i70 = 0; i70 < nlmmt; i70++) wk[i70] = w[i70][1] * fam;
	// returns
	skip62 = true;
      }
    }
    if (!skip62) {
      if (lmode == 0 || lmode == 1) { // label 62
	for (int i63 = 0; i63 < nlmmt; i63++) wk[i63] = w[i63][0] * cf1 + w[i63][1] * cf2;
      }
    }
  }
  // Clean up memory
  delete[] up;
  delete[] un;
  for (int wi = nlmmt - 1; wi > -1; wi--) delete[] w[wi];
  delete[] w;
  delete[] ylm;
}
