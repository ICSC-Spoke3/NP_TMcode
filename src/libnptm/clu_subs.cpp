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

/*! \file clu_subs.cpp
 *
 * \brief C++ implementation of CLUSTER subroutines.
 */

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

#ifndef INCLUDE_CLU_SUBS_H_
#include "../include/clu_subs.h"
#endif

#ifdef USE_NVTX
#include <nvtx3/nvToolsExt.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void apc(
	 double ****zpv, int le, dcomplex **am0m, dcomplex **w,
	 double sqk, double **gapr, dcomplex **gapp
) {
  dcomplex **ac, **gap;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex uimmp, summ, sume, suem, suee, summp, sumep;
  dcomplex suemp, sueep;
  double cof = 1.0 / sqk;
  double cimu = cof / sqrt(2.0);
  int nlem = le * (le + 2);
  const int nlemt = nlem + nlem;
  ac = new dcomplex*[nlemt];
  gap = new dcomplex*[3];
  for (int ai = 0; ai < nlemt; ai++) ac[ai] = new dcomplex[2]();
  for (int gi = 0; gi < 3; gi++) gap[gi] = new dcomplex[2]();
  for (int j45 = 1; j45 <= nlemt; j45++) {
    int j = j45 - 1;
    ac[j][0] = cc0;
    ac[j][1] = cc0;
    for (int i45 = 1; i45 <= nlemt; i45++) {
      int i = i45 - 1;
      ac[j][0] += (am0m[j][i] * w[i][0]);
      ac[j][1] += (am0m[j][i] * w[i][1]);
    } //i45 loop
  } //j45 loop
  for (int imu90 = 1; imu90 <=3; imu90++) {
    int mu = imu90 - 2;
    gap[imu90 - 1][0] = cc0;
    gap[imu90 - 1][1] = cc0;
    gapp[imu90 - 1][0] = cc0;
    gapp[imu90 - 1][1] = cc0;
    for (int l80 =1; l80 <= le; l80++) {
      int lpo = l80 + 1;
      int ltpo = lpo + l80;
      int imm = l80 * lpo;
      for (int ilmp = 1; ilmp <= 3; ilmp++) {
	if ((l80 == 1 && ilmp == 1) || (l80 == le && ilmp == 3)) continue; // ilmp loop
	int lmpml = ilmp - 2;
	int lmp = l80 + lmpml;
	uimmp = (-1.0 * lmpml) * uim;
	int impmmmp = lmp * (lmp + 1);
	for (int im70 = 1; im70 <= ltpo; im70++) {
	  int m = im70 - lpo;
	  int mmp = m - mu;
	  int abs_mmp = (mmp > 0) ? mmp : -mmp;
	  if (abs_mmp <= lmp) {
	    int i = imm + m;
	    int ie = i + nlem;
	    int imp = impmmmp + mmp;
	    int impe = imp + nlem;
	    double cgc = cg1(lmpml, mu, l80, m);
	    int jpo = 2;
	    for (int ipo = 1; ipo <= 2; ipo++) {
	      if (ipo == 2) jpo = 1;
	      summ = dconjg(ac[i - 1][ipo - 1]) * ac[imp - 1][ipo - 1];
	      sume = dconjg(ac[i - 1][ipo - 1]) * ac[impe - 1][ipo - 1];
	      suem = dconjg(ac[ie - 1][ipo - 1]) * ac[imp - 1][ipo - 1];
	      suee = dconjg(ac[ie - 1][ipo - 1]) * ac[impe - 1][ipo - 1];
	      summp = dconjg(ac[i - 1][jpo - 1]) * ac[imp - 1][ipo - 1];
	      sumep = dconjg(ac[i - 1][jpo - 1]) * ac[impe - 1][ipo - 1];
	      suemp = dconjg(ac[ie - 1][jpo - 1]) * ac[imp - 1][ipo - 1];
	      sueep = dconjg(ac[ie - 1][jpo - 1]) * ac[impe - 1][ipo - 1];
	      if (lmpml != 0) {
		summ *= uimmp;
		sume *= uimmp;
		suem *= uimmp;
		suee *= uimmp;
		summp *= uimmp;
		sumep *= uimmp;
		suemp *= uimmp;
		sueep *= uimmp;
	      }
	      // label 55
	      gap[imu90 - 1][ipo - 1] += (
					  (
					   summ * zpv[l80 - 1][ilmp - 1][0][0]
					   + sume * zpv[l80 - 1][ilmp - 1][0][1]
					   + suem * zpv[l80 - 1][ilmp - 1][1][0]
					   + suee * zpv[l80 - 1][ilmp - 1][1][1]
					   ) * cgc
					  );
	      gapp[imu90 - 1][ipo - 1] += (
					   (
					    summp * zpv[l80 - 1][ilmp - 1][0][0]
					    + sumep * zpv[l80 - 1][ilmp - 1][0][1]
					    + suemp * zpv[l80 - 1][ilmp - 1][1][0]
					    + sueep * zpv[l80 - 1][ilmp - 1][1][1]
					    ) * cgc
					   );
	    } // ipo loop
	  } // ends im70 loop
	} // im70 loop
      } // ilmp loop
    } // l80 loop
  } // imu90 loop
  for (int ipo95 = 1; ipo95 <= 2; ipo95++) {
    sume = gap[0][ipo95 - 1] * cimu;
    suee = gap[1][ipo95 - 1] * cof;
    suem = gap[2][ipo95 - 1] * cimu;
    gapr[0][ipo95 - 1] = real(sume - suem);
    gapr[1][ipo95 - 1] = real((sume + suem) * uim);
    gapr[2][ipo95 - 1] = real(suee);
    sumep = gapp[0][ipo95 - 1] * cimu;
    sueep = gapp[1][ipo95 - 1] * cof;
    suemp = gapp[2][ipo95 - 1] * cimu;
    gapp[0][ipo95 - 1] = sumep - suemp;
    gapp[1][ipo95 - 1] = (sumep + suemp) * uim;
    gapp[2][ipo95 - 1] = sueep;
  } // ipo95 loop
  // Clean memory
  for (int ai = nlemt - 1; ai > -1; ai--) delete[] ac[ai];
  for (int gi = 2; gi > -1; gi--) delete[] gap[gi];
  delete[] ac;
  delete[] gap;
}

void apcra(
	   double ****zpv, const int le, dcomplex **am0m, int inpol, double sqk,
	   double **gaprm, dcomplex **gappm
) {
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex uimtl, uimtls, ca11, ca12, ca21, ca22;
  dcomplex a11, a12, a21, a22, sum1, sum2, fc;
  double ****svw = new double***[le];
  dcomplex ****svs = new dcomplex***[le];
  for (int i = 0; i < le; i++) {
    svw[i] = new double**[3];
    svs[i] = new dcomplex**[3];
    for (int j = 0; j < 3; j++) {
      svw[i][j] = new double*[2];
      svs[i][j] = new dcomplex*[2];
      for (int k = 0; k < 2; k++) {
	svw[i][j][k] = new double[2]();
	svs[i][j][k] = new dcomplex[2]();
      }
    }
  }
  int nlem = le * (le + 2);
  for (int l28 = 1; l28 <= le; l28++) {
    int lpo = l28 + 1;
    int ltpo = lpo + l28;
    double fl = sqrt(1.0 * ltpo);
    for (int ilmp = 1; ilmp <= 3; ilmp++) {
      if ((l28 == 1 && ilmp == 1) || (l28 == le && ilmp == 3)) continue; // ilmp loop
      int lmpml = ilmp - 2;
      int lmp = l28 + lmpml;
      double flmp = sqrt(1.0 * (lmp + lmp + 1));
      double fllmp = flmp / fl;
      double cgmmo = fllmp * cg1(lmpml, 0, l28, 1);
      double cgmpo = fllmp * cg1(lmpml, 0, l28, -1);
      if (inpol == 0) {
	double cgs = cgmpo + cgmmo;
	double cgd = cgmpo - cgmmo;
	svw[l28 - 1][ilmp - 1][0][0] = cgs;
	svw[l28 - 1][ilmp - 1][0][1] = cgd;
	svw[l28 - 1][ilmp - 1][1][0] = cgd;
	svw[l28 - 1][ilmp - 1][1][1] = cgs;
      } else { // label 22
	svw[l28 - 1][ilmp - 1][0][0] = cgmpo;
	svw[l28 - 1][ilmp - 1][1][0] = cgmpo;
	svw[l28 - 1][ilmp - 1][0][1] = -cgmmo;
	svw[l28 - 1][ilmp - 1][1][1] = cgmmo;
      }
      // label 26
    } // ilmp loop
  } // l28 loop
  for (int l30 = 1; l30 <= le; l30++) { // 0-init: can be omitted
    for (int ilmp = 1; ilmp <= 3; ilmp++) {
      for (int ipa = 1; ipa <= 2; ipa++) {
	for (int ipamp = 1; ipamp <= 2; ipamp++) {
	  svs[l30 - 1][ilmp - 1][ipa - 1][ipamp - 1] = cc0;
	}
      } // ipa loop
    } // ilmp loop
  } // l30 loop
  for (int l58 = 1; l58 <= le; l58 ++) {
    int lpo = l58 + 1;
    int ltpo = l58 + lpo;
    int imm = l58 * lpo;
    for (int ilmp = 1; ilmp <= 3; ilmp++) {
      if ((l58 == 1 && ilmp == 1) || (l58 == le && ilmp == 3)) continue; // ilmp loop
      int lmpml = ilmp - 2;
      int lmp = l58 + lmpml;
      int impmm = lmp * (lmp + 1);
      uimtl = uim * (1.0 * lmpml);
      if (lmpml == 0) uimtl = 1.0 + 0.0 * I;
      for (int im54 = 1; im54 <= ltpo; im54++) {
	int m = im54 - lpo;
	int i = imm + m;
	int ie = i + nlem;
	for (int imu52 = 1; imu52 <= 3; imu52++) {
	  int mu = imu52 - 2;
	  int mmp = m - mu;
	  int abs_mmp = (mmp > 0) ? mmp : -mmp;
	  if (abs_mmp <= lmp) {
	    int imp = impmm + mmp;
	    int impe = imp + nlem;
	    double cgc = cg1(lmpml, -mu, l58, -m);
	    for (int ls = 1; ls <= le; ls++) {
	      int lspo = ls + 1;
	      int lstpo = ls + lspo;
	      int ismm = ls * lspo;
	      for (int ilsmp = 1; ilsmp <= 3; ilsmp++) {
		if ((ls == 1 && ilsmp == 1) || (ls == le && ilsmp == 3)) continue; // ilsmp loop
		int lsmpml = ilsmp - 2;
		int lsmp = ls + lsmpml;
		int ismpmm = lsmp * (lsmp + 1);
		uimtls = -uim * (1.0 * lsmpml);
		if (lsmpml == 0) uimtls = 1.0 + 0.0 * I;
		for (int ims = 1; ims <= lstpo; ims++) {
		  int ms = ims - lspo;
		  int msmp = ms - mu;
		  int abs_msmp = (msmp > 0) ? msmp : -msmp;
		  if (abs_msmp <= lsmp) {
		    int is = ismm + ms;
		    int ise = is + nlem;
		    int ismp = ismpmm + msmp;
		    int ismpe = ismp + nlem;
		    double cgcs = cg1(lsmpml, mu, ls, ms);
		    fc = (uimtl * uimtls) * (cgc * cgcs);
		    ca11 = dconjg(am0m[is - 1][i - 1]);
		    ca12 = dconjg(am0m[is - 1][ie - 1]);
		    ca21 = dconjg(am0m[ise - 1][i - 1]);
		    ca22 = dconjg(am0m[ise - 1][ie - 1]);
		    a11 = am0m[ismp - 1][imp - 1];
		    a12 = am0m[ismp - 1][impe - 1];
		    a21 = am0m[ismpe - 1][imp - 1];
		    a22 = am0m[ismpe - 1][impe - 1];
		    double z11 = zpv[ls - 1][ilsmp - 1][0][0];
		    double z12 = zpv[ls - 1][ilsmp - 1][0][1];
		    double z21 = zpv[ls - 1][ilsmp - 1][1][0];
		    double z22 = zpv[ls - 1][ilsmp - 1][1][1];
		    svs[l58 - 1][ilmp - 1][0][0] += ((ca11 * a11 * z11
						      + ca11 * a21 * z12
						      + ca21 * a11 * z21
						      + ca21 * a21 * z22) * fc);
		    svs[l58 - 1][ilmp - 1][0][1] += ((ca11 * a12 * z11
						      + ca11 * a22 * z12
						      + ca21 * a12 * z21
						      + ca21 * a22 * z22) * fc);
		    svs[l58 - 1][ilmp - 1][1][0] += ((ca12 * a11 * z11
						      + ca12 * a21 * z12
						      + ca22 * a11 * z21
						      + ca22 * a21 * z22) * fc);
		    svs[l58 - 1][ilmp - 1][1][1] += ((ca12 * a12 * z11
						      + ca12 * a22 * z12
						      + ca22 * a12 * z21
						      + ca22 * a22 * z22) * fc);
		  } // ends ims loop
		} // ims loop
	      } // ilsmp loop
	    } // ls loop
	  } // ends imu52 loop
	} // imu52 loop
      } // im54 loop
    } // ilmp loop
  } // l58 loop
  sum1 = cc0;
  sum2 = cc0;
  for (int l68 = 1; l68 <= le; l68++) {
    //int lpo = l68 + 1;
    //int ltpo = l68 + lpo;
    //int imm = l68 * lpo;
    for (int ilmp = 1; ilmp <= 3; ilmp++) {
      if ((l68 == 1 && ilmp == 1) || (l68 == le && ilmp == 3)) continue; // ilmp loop
      if (inpol == 0) {
	sum1 += (
		 svw[l68 - 1][ilmp - 1][0][0] * svs[l68 - 1][ilmp - 1][0][0]
		 + svw[l68 - 1][ilmp - 1][1][0] * svs[l68 - 1][ilmp - 1][0][1]
		 + svw[l68 - 1][ilmp - 1][1][0] * svs[l68 - 1][ilmp - 1][1][0]
		 + svw[l68 - 1][ilmp - 1][0][0] * svs[l68 - 1][ilmp - 1][1][1]
		 );
	sum2 += (
		 svw[l68 - 1][ilmp - 1][0][1] * svs[l68 - 1][ilmp - 1][0][0]
		 + svw[l68 - 1][ilmp - 1][1][1] * svs[l68 - 1][ilmp - 1][0][1]
		 + svw[l68 - 1][ilmp - 1][1][1] * svs[l68 - 1][ilmp - 1][1][0]
		 + svw[l68 - 1][ilmp - 1][0][1] * svs[l68 - 1][ilmp - 1][1][1]
		 );
      } else { // label 62
	sum1 += (
		 svw[l68 - 1][ilmp - 1][1][0] * svs[l68 - 1][ilmp - 1][0][0]
		 + svw[l68 - 1][ilmp - 1][0][0] * svs[l68 - 1][ilmp - 1][0][1]
		 + svw[l68 - 1][ilmp - 1][0][0] * svs[l68 - 1][ilmp - 1][1][0]
		 + svw[l68 - 1][ilmp - 1][1][0] * svs[l68 - 1][ilmp - 1][1][1]
		 );
	sum2 += (
		 svw[l68 - 1][ilmp - 1][1][1] * svs[l68 - 1][ilmp - 1][0][0]
		 + svw[l68 - 1][ilmp - 1][0][1] * svs[l68 - 1][ilmp - 1][0][1]
		 + svw[l68 - 1][ilmp - 1][0][1] * svs[l68 - 1][ilmp - 1][1][0]
		 + svw[l68 - 1][ilmp - 1][1][1] * svs[l68 - 1][ilmp - 1][1][1]
		 );
      } // label 66, ends ilmp loop
    } // ilmp loop
  } // l68 loop
  const double half_pi = acos(0.0);
  double cofs = half_pi * 2.0 / sqk;
  gaprm[0][0] = 0.0;
  gaprm[0][1] = 0.0;
  gaprm[1][0] = 0.0;
  gaprm[1][1] = 0.0;
  gappm[0][0] = cc0;
  gappm[0][1] = cc0;
  gappm[1][0] = cc0;
  gappm[1][1] = cc0;
  if (inpol == 0) {
    sum1 *= cofs;
    sum2 *= cofs;
    gaprm[2][0] = real(sum1);
    gaprm[2][1] = real(sum1);
    gappm[2][0] = sum2 * uim;
    gappm[2][1] = -gappm[2][0];
  } else { // label 72
    cofs *= 2.0;
    gaprm[2][0] = real(sum1) * cofs;
    gaprm[2][1] = real(sum2) * cofs;
    gappm[2][0] = cc0;
    gappm[2][1] = cc0;
  }
  
  // Clean memory
  for (int i = le - 1; i > -1; i--) {
    for (int j = 2; j > -1; j--) {
      for (int k = 1; k > -1; k--) {
	delete[] svw[i][j][k];
	delete[] svs[i][j][k];
      }
      delete[] svw[i][j];
      delete[] svs[i][j];
    }
    delete[] svw[i];
    delete[] svs[i];
  }
  delete[] svw;
  delete[] svs;
}

dcomplex cdtp(dcomplex z, dcomplex **am, int i, int jf, int k, int nj) {
  /* NOTE: the original FORTRAN code treats the AM matrix as a
   * vector. This is not directly allowed in C++ and it requires
   * accounting for the different dimensions.
   */
  dcomplex result = z;
  if (nj > 0) {
    int jl = jf + nj - 1;
    for (int j = jf; j <= jl; j++) {
      result += (am[i - 1][j - 1] * am[j - 1][k - 1]);
    }
  }
  return result;
}

// #ifdef USE_TARGET_OFFLOAD
// #pragma omp begin declare target device_type(any)
// #endif
double cgev(int ipamo, int mu, int l, int m) {
  double result = 0.0;
  double xd = 0.0, xn = 0.0;
  if (ipamo == 0) {
    if (m != 0 || mu != 0) { // label 10
      if (mu != 0) {
	xd = 2.0 * l * (l + 1);
	if (mu <= 0) {
	  xn = 1.0 * (l + m) * (l - m + 1);
	  result = sqrt(xn / xd);
	} else { // label 15
	  xn = 1.0 * (l - m) * (l + m + 1);
	  result = -sqrt(xn / xd);
	}
      } else { // label 20
	xd = 1.0 * (l + 1) * l;
	xn = -1.0 * m;
	result = xn / sqrt(xd);
      }
    }
  } else { // label 30
    xd = 2.0 * l * (l * 2 - 1);
    if (mu < 0) { // label 35
      xn = 1.0 * (l - 1 + m) * (l + m);
    } else if (mu == 0) { // label 40
      xn = 2.0 * (l - m) * (l + m);
    } else { // mu > 0, label 45
      xn = 1.0 * (l - 1 - m) * (l - m);
    }
    result = sqrt(xn / xd);
  }
  return result;
}
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp end declare target
// #endif

void cms(dcomplex **am, ParticleDescriptor *c1) {
  dcomplex dm, de, cgh, cgk;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  int ndi = c1->nsph * c1->nlim;
  int nbl = 0;
  int nsphmo = c1->nsph - 1;
  for (int n1 = 1; n1 <= nsphmo; n1++) { // GPU portable?
    int in1 = (n1 - 1) * c1->nlim;
    int n1po = n1 + 1;
    for (int n2 = n1po; n2 <= c1->nsph; n2++) {
      int in2 = (n2 - 1) * c1->nlim;
      nbl++;
      for (int l1 = 1; l1 <= c1->li; l1++) {
	int l1po = l1 + 1;
	int il1 = l1po * l1;
	int l1tpo = l1po + l1;
	for (int im1 = 1; im1 <= l1tpo; im1++) {
	  int m1 = im1 - l1po;
	  int ilm1 = il1 + m1;
	  int ilm1e = ilm1 + ndi;
	  int i1 = in1 + ilm1;
	  int i1e = in1 + ilm1e;
	  int j1 = in2 + ilm1;
	  int j1e = in2 + ilm1e;
	  for (int l2 = 1; l2 <= c1->li; l2++) {
	    int l2po = l2 + 1;
	    int il2 = l2po * l2;
	    int l2tpo = l2po + l2;
	    int ish = ((l2 + l1) % 2 == 0) ? 1 : -1;
	    int isk = -ish;
	    for (int im2 = 1; im2 <= l2tpo; im2++) {
	      int m2 = im2 - l2po;
	      int ilm2 = il2 + m2;
	      int ilm2e = ilm2 + ndi;
	      int i2 = in2 + ilm2;
	      int i2e = in2 + ilm2e;
	      int j2 = in1 + ilm2;
	      int j2e = in1 + ilm2e;
	      cgh = ghit(0, 0, nbl, l1, m1, l2, m2, c1);
	      cgk = ghit(0, 1, nbl, l1, m1, l2, m2, c1);
	      am[i1 - 1][i2 - 1] = cgh;
	      am[i1 - 1][i2e - 1] = cgk;
	      am[i1e - 1][i2 - 1] = cgk;
	      am[i1e - 1][i2e - 1] = cgh;
	      am[j1 - 1][j2 - 1] = cgh * (1.0 * ish);
	      am[j1 - 1][j2e - 1] = cgk * (1.0 * isk);
	      am[j1e - 1][j2 - 1] = cgk * (1.0 * isk);
	      am[j1e - 1][j2e - 1] = cgh * (1.0 * ish);
	    }
	  }
	} // im1 loop
      } // l1 loop
    } // n2 loop
  } // n1 loop
  for (int n1 = 1; n1 <= c1->nsph; n1++) { // GPU portable?
    int in1 = (n1 - 1) * c1->nlim;
    for (int l1 = 1; l1 <= c1->li; l1++) {
      dm = c1->rmi[l1 - 1][n1 - 1];
      de = c1->rei[l1 - 1][n1 - 1];
      int l1po = l1 + 1;
      int il1 = l1po * l1;
      int l1tpo = l1po + l1;
      for (int im1 = 1; im1 <= l1tpo; im1++) {
	int m1 = im1 - l1po;
	int ilm1 = il1 + m1;
	int i1 = in1 + ilm1;
	int i1e = i1 + ndi;
	for (int ilm2 = 1; ilm2 <= c1->nlim; ilm2++) {
	  int i2 = in1 + ilm2;
	  int i2e = i2 + ndi;
	  am[i1 - 1][i2 - 1] = cc0;
	  am[i1 - 1][i2e - 1] = cc0;
	  am[i1e - 1][i2 - 1] = cc0;
	  am[i1e - 1][i2e - 1] = cc0;
	}
	am[i1 - 1][i1 - 1] = dm;
	am[i1e - 1][i1e - 1] = de;
      } // im1 loop
    } // l1 loop
  } // n1 loop
}

void crsm1(double vk, double exri, ParticleDescriptor *c1) {
  dcomplex ***svf, ***svw, **svs;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  dcomplex cam = cc0;
  const int le4po = 4 * c1->le + 1;
  svf = new dcomplex**[le4po];
  svw = new dcomplex**[le4po];
  svs = new dcomplex*[le4po];
  for (int si = 0; si < le4po; si++) {
    svf[si] = new dcomplex*[le4po];
    svw[si] = new dcomplex*[4];
    svs[si] = new dcomplex[4]();
    for (int sj = 0; sj < le4po; sj++) svf[si][sj] = new dcomplex[4]();
    for (int sj = 0; sj < 4; sj++) svw[si][sj] = new dcomplex[4]();
  }
  double exdc = exri * exri;
  double ccs = 1.0 / (vk * vk);
  const double pi4sq = 64.0 * acos(0.0) * acos(0.0);
  double cint = ccs / (pi4sq * exdc);
  int letpo = c1->le + c1->le + 1;
  for (int i20 = 0; i20 < 16; i20++) c1->vintm[i20] = cc0; // 0-init: can be omitted
  for (int lpo40 = 1; lpo40 <= letpo; lpo40++) {
    int l = lpo40 - 1;
    int ltpo = lpo40 + l;
    int immn = letpo - l;
    int immx = letpo + l;
    for (int imf = immn; imf <= immx; imf++) { // 0-init: can be omitted
      for (int ims = immn; ims <= immx; ims++) {
	for (int ipo = 1; ipo <= 4; ipo++) {
	  svf[imf - 1][ims - 1][ipo - 1] = cc0;
	} // ipo loop
      } // ims loop
    } // imf loop
    for (int l1 = 1; l1 <= c1->le; l1++) {
      int il1 = l1 * (l1 + 1);
      for (int l2 = 1; l2 <= c1->le; l2++) {
	int abs_l2ml1 = (l2 > l1) ? l2 - l1 : l1 - l2;
	if (l < abs_l2ml1 || l > l2 + l1) continue; // l2 loop
	int il2 = l2 * (l2 + 1);
	for (int im = immn; im >= immx; im++) { // 0-init: can be omitted
	  for (int ipa = 1; ipa <= 4; ipa++) {
	    svs[im - 1][ipa - 1] = cc0;
	    for (int ipo = 1; ipo <= 4; ipo++) {
	      svw[im - 1][ipa - 1][ipo - 1] = cc0;
	    } // ipo loop
	  } // ipa loop
	} // im loop
	for (int im = immn; im <= immx; im++) {
	  int m = im - letpo;
	  r3jmr(l, l1, l2, m, c1->rac3j);
	  int m1mnmo = (-l1 > -l2 - m) ? -(l1 + 1) : -(l2 + m + 1);
	  int nm1 = (l1 < l2 - m) ? (l1 - m1mnmo) : (l2 - m - m1mnmo);
	  for (int im1 = 1; im1 <= nm1; im1++) {
	    int m1 = -im1 - m1mnmo;
	    int isn = 1;
	    if (m1 % 2 != 0) isn = -1;
	    double cg3j = c1->rac3j[im1 - 1] * isn;
	    int ilm1 = il1 + m1;
	    int ilm2 = il2 + m1 - m;
	    int ipa = 0;
	    for (int ipa1 = 1; ipa1 <= 2; ipa1++) {
	      int i1 = ilm1;
	      if (ipa1 == 2) i1 = ilm1 + c1->nlem;
	      for (int ipa2 = 1; ipa2 <= 2; ipa2++) {
		int i2 = ilm2;
		if (ipa2 == 2) i2 = ilm2 + c1->nlem;
		ipa++;
		svs[im - 1][ipa - 1] += (c1->am0m[i1 - 1][i2 - 1] * cg3j);
		int ipo = 0;
		for (int ipo2 = 1; ipo2 <= 2; ipo2++) {
		  for (int ipo1 = 3; ipo1 <= 4; ipo1++) {
		    ipo++;
		    svw[im - 1][ipa - 1][ipo - 1] += (c1->w[i1 - 1][ipo1 - 1] * c1->w[i2 - 1][ipo2 - 1] * cg3j);
		  } // ipo1 loop
		} // ipo2 loop
	      } // ipa2 loop
	    } // ipa1 loop
	  } // im1 loop
	  // label 32 loops
	  for (int imf = immn; imf <= immx; imf++) {
	    for (int ims = immn; ims <= immx; ims++) {
	      for (int ipo = 1; ipo <= 4; ipo++) {
		for (int ipa = 1; ipa <= 4; ipa++) {
		  svf[imf - 1][ims - 1][ipo - 1] += (svw[imf - 1][ipa - 1][ipo - 1] * svs[ims - 1][ipa - 1]);
		} // ipa loop
	      } // ipo loop
	    } // ims loop
	  } // imf loop
	  // ends loop level 34, which are l2 loop and l1 loop
	} // im loop
      } // l2 loop
    } // l1 loop
    for (int imf = immn; imf <= immx; imf++) {
      for (int ims = immn; ims <= immx; ims++) {
	int i = 0;
	for (int ipo1 = 1; ipo1 <= 4; ipo1++) {
	  cam = dconjg(svf[imf - 1][ims - 1][ipo1 - 1]);
	  for (int ipo2 = 1; ipo2 <= 4; ipo2++) {
	    i++;
	    c1->vintm[i - 1] += (svf[imf - 1][ims - 1][ipo2 - 1] * cam * (1.0 * ltpo));
	  }
	} // ipo1 loop
      } // ims loop
    } // imf loop
  } // lpo40 loop
  for (int i42 = 0; i42 < 16; i42++) c1->vintm[i42] *= cint;

  // Clean memory
  for (int si = le4po - 1; si > -1; si--) {
    for (int sj = le4po - 1; sj > -1; sj--) delete[] svf[si][sj];
    for (int sj = 3; sj > -1; sj--) delete[] svw[si][sj];
    delete[] svf[si];
    delete[] svw[si];
    delete[] svs[si];
  }
  delete[] svf;
  delete[] svw;
  delete[] svs;
}

// #ifdef USE_TARGET_OFFLOAD
// #pragma omp begin declare target device_type(any)
// #endif
dcomplex ghit_d(
	      int ihi, int ipamo, int nbl, int l1, int m1, int l2, int m2,
	      ParticleDescriptor *c1, double *rac3j
) {
  /* NBL identifies transfer vector going from N2 to N1;
   * IHI=0 for Hankel, IHI=1 for Bessel, IHI=2 for Bessel from origin;
   * depending on IHI, IPAM=0 gives H or I, IPAM= 1 gives K or L. */
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex csum = cc0, cfun = cc0;
  dcomplex result = cc0;

  if (ihi == 2) {
    if (c1->rxx[nbl - 1] == 0.0 && c1->ryy[nbl - 1] == 0.0 && c1->rzz[nbl - 1] == 0.0) {
      if (ipamo == 0) {
	if (l1 == l2 && m1 == m2) result = 1.0 + 0.0 * I;
      }
      return result;
    }
  }
  // label 10
  int l1mp = l1 - ipamo;
  int l1po = l1 + 1;
  int m1mm2 = m1 - m2;
  int m1mm2m = (m1mm2 > 0) ? m1mm2 + 1 : 1 - m1mm2;
  int lminpo = (l2 - l1mp > 0) ? l2 - l1mp + 1 : l1mp - l2 + 1;
  int lmaxpo = l2 + l1mp + 1;
  int i3j0in = c1->ind3j[l1mp][l2 - 1];
  int ilin = -1;
  if (m1mm2m > lminpo && (m1mm2m - lminpo) % 2 != 0) ilin = 0;
  int isn = 1;
  if (m1 % 2 != 0) isn *= -1;
  if (lminpo % 2 == 0) {
    isn *= -1;
    if (l2 > l1mp) isn *= -1;
  }
  // label 12
  int nblmo = nbl - 1;
  if (ihi != 2) {
    int nbhj = nblmo * c1->litpo;
    int nby = nblmo * c1->litpos;
    if (ihi != 1) {
      for (int jm24 = 1; jm24 <= 3; jm24++) {
	csum = cc0;
	int mu = jm24 - 2;
	int mupm1 = mu + m1;
	int mupm2 = mu + m2;
	if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	  int jsn = -isn;
	  if (mu == 0) jsn = isn;
	  double cr = cgev(ipamo, mu, l1, m1) * cgev(0, mu, l2, m2);
	  int i3j0 = i3j0in;
	  if (mupm1 == 0 && mupm2 == 0) {
	    int lt14 = lminpo;
	    while (lt14 <= lmaxpo) {
	      i3j0++;
	      int l3 = lt14 - 1;
	      int ny = l3 * l3 + lt14;
	      double aors = 1.0 * (l3 + lt14);
	      double f3j = (c1->v3j0[i3j0 - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	      cfun = (c1->vh[nbhj + lt14 - 1] * c1->vyhj[nby + ny - 1]) * f3j;
	      csum += cfun;
	      jsn *= -1;
	      lt14 += 2;
	    }
	    // goes to 22
	  } else { // label 16
	    r3jjr_d(l1mp, l2, -mupm1, mupm2, rac3j);
	    int il = ilin;
	    int lt20 = lminpo;
	    while (lt20 <= lmaxpo) {
	      i3j0++;
	      if (m1mm2m <= lt20) {
		il += 2;
		int l3 = lt20 - 1;
		int ny = l3 * l3  + lt20 + m1mm2;
		double aors = 1.0 * (l3 + lt20);
		double f3j = (rac3j[il - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
		cfun = (c1->vh[nbhj + lt20 - 1] * c1->vyhj[nby + ny - 1]) * f3j;
		csum += cfun;
	      }
	      // label 20
	      jsn *= -1;
	      lt20 += 2;
	    }
	  }
	  // label 22
	  csum *= cr;
	  result += csum;
	}
	// Otherwise there is nothing to add
      } // jm24 loop. Should go to 70
    } else { // label 30, IHI == 1
      for (int jm44 = 1; jm44 <= 3; jm44++) {
	csum = cc0;
	int mu = jm44 - 2;
	int mupm1 = mu + m1;
	int mupm2 = mu + m2;
	if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	  int jsn = - isn;
	  if (mu == 0) jsn = isn;
	  double cr = cgev(ipamo, mu, l1, m1) * cgev(0, mu, l2, m2);
	  int i3j0 = i3j0in;
	  if (mupm1 == 0 && mupm2 == 0) {
	    int lt34 = lminpo;
	    while (lt34 <= lmaxpo) {
	      i3j0++;
	      int l3 = lt34 - 1;
	      int ny = l3 * l3 + lt34;
	      double aors = 1.0 * (l3 + lt34);
	      double f3j = (c1->v3j0[i3j0 - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	      cfun = (c1->vh[nbhj + lt34 - 1] * c1->vyhj[nby + ny - 1]) * f3j;
	      csum += cfun;
	      jsn *= -1;
	      lt34 += 2;
	    }
	    // goes to 42
	  } else { // label 36
	    r3jjr_d(l1mp, l2, -mupm1, mupm2, rac3j);
	    int il = ilin;
	    int lt40 = lminpo;
	    while (lt40 <= lmaxpo) {
	      i3j0++;
	      if (m1mm2m <= lt40) {
		il += 2;
		int l3 = lt40 - 1;
		int ny = l3 * l3  + lt40 + m1mm2;
		double aors = 1.0 * (l3 + lt40);
		double f3j = (rac3j[il - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
		cfun = (c1->vh[nbhj + lt40 - 1] * c1->vyhj[nby + ny - 1]) * f3j;
		csum += cfun;
	      }
	      // label 40
	      jsn *= -1;
	      lt40 += 2;
	    }
	  }
	  // label 42
	  csum *= cr;
	  result += csum;
	}
	// Otherwise there is nothing to add
      } // jm44 loop. Should go to 70
    }
    // goes to 70
  } else { // label 50, IHI == 2
    int nbhj = nblmo * c1->lmtpo;
    int nby = nblmo * c1->lmtpos;
    for (int jm64 = 1; jm64 <= 3; jm64++) {
      csum = cc0;
      int mu = jm64 - 2;
      int mupm1 = mu + m1;
      int mupm2 = mu + m2;
      if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	int jsn = -isn;
	if (mu == 0) jsn = isn;
	double cr = cgev(ipamo, mu, l1, m1) * cgev(0, mu, l2, m2);
	int i3j0 = i3j0in;
	if (mupm1 == 0 && mupm2 == 0) {
	  int lt54 = lminpo;
	  while (lt54 <= lmaxpo) {
	    i3j0++;
	    int l3 = lt54 - 1;
	    int ny = l3 * l3 + lt54;
	    double aors = 1.0 * (l3 + lt54);
	    double f3j = (c1->v3j0[i3j0 - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	    cfun = (c1->vj0[nbhj + lt54 - 1] * c1->vyj0[nby + ny - 1]) * f3j;
	    csum += cfun;
	    jsn *= -1;
	    lt54 += 2;
	  }
	  // goes to 62
	} else { // label 56
	  r3jjr_d(l1mp, l2, -mupm1, mupm2, rac3j);
	  int il = ilin;
	  int lt60 = lminpo;
	  while (lt60 <= lmaxpo) {
	    i3j0++;
	    if (m1mm2m <= lt60) {
	      il += 2;
	      int l3 = lt60 - 1;
	      int ny = l3 * l3  + lt60 + m1mm2;
	      double aors = 1.0 * (l3 + lt60);
	      double f3j = (rac3j[il - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	      cfun = (c1->vj0[nbhj + lt60 - 1] * c1->vyj0[nby + ny - 1]) * f3j;
	      csum += cfun;
	    }
	    // label 60
	    jsn *= -1;
	    lt60 += 2;
	  }
	}
	// label 62
	csum *= cr;
	result += csum;
      }
      // Otherwise there is nothing to add
    } // jm64 loop. Should go to 70
  }
  // label 70
  const double four_pi = acos(0.0) * 8.0;
  if (ipamo != 1) {
    double cr = sqrt(four_pi * (l1 + l1po) * (l2 + l2 + 1));
    result *= cr;
  } else {
    double cr = sqrt(four_pi * (l1 + l1mp) * (l1 + l1po) * (l2 + l2 + 1) / l1po);
    result *= (cr * uim);
  }
  return result;
}
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp end declare target
// #endif

// #ifdef USE_TARGET_OFFLOAD
// #pragma omp begin declare target device_type(any)
// #endif
dcomplex ghit(
	      int ihi, int ipamo, int nbl, int l1, int m1, int l2, int m2,
	      ParticleDescriptor *c1
) {
  /* NBL identifies transfer vector going from N2 to N1;
   * IHI=0 for Hankel, IHI=1 for Bessel, IHI=2 for Bessel from origin;
   * depending on IHI, IPAM=0 gives H or I, IPAM= 1 gives K or L. */
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex csum = cc0, cfun = cc0;
  dcomplex result = cc0;

  if (ihi == 2) {
    if (c1->rxx[nbl - 1] == 0.0 && c1->ryy[nbl - 1] == 0.0 && c1->rzz[nbl - 1] == 0.0) {
      if (ipamo == 0) {
	if (l1 == l2 && m1 == m2) result = 1.0 + 0.0 * I;
      }
      return result;
    }
  }
  // label 10
  int l1mp = l1 - ipamo;
  int l1po = l1 + 1;
  int m1mm2 = m1 - m2;
  int m1mm2m = (m1mm2 > 0) ? m1mm2 + 1 : 1 - m1mm2;
  int lminpo = (l2 - l1mp > 0) ? l2 - l1mp + 1 : l1mp - l2 + 1;
  int lmaxpo = l2 + l1mp + 1;
  int i3j0in = c1->ind3j[l1mp][l2 - 1];
  int ilin = -1;
  if (m1mm2m > lminpo && (m1mm2m - lminpo) % 2 != 0) ilin = 0;
  int isn = 1;
  if (m1 % 2 != 0) isn *= -1;
  if (lminpo % 2 == 0) {
    isn *= -1;
    if (l2 > l1mp) isn *= -1;
  }
  // label 12
  int nblmo = nbl - 1;
  if (ihi != 2) {
    int nbhj = nblmo * c1->litpo;
    int nby = nblmo * c1->litpos;
    if (ihi != 1) {
      for (int jm24 = 1; jm24 <= 3; jm24++) {
	csum = cc0;
	int mu = jm24 - 2;
	int mupm1 = mu + m1;
	int mupm2 = mu + m2;
	if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	  int jsn = -isn;
	  if (mu == 0) jsn = isn;
	  double cr = cgev(ipamo, mu, l1, m1) * cgev(0, mu, l2, m2);
	  int i3j0 = i3j0in;
	  if (mupm1 == 0 && mupm2 == 0) {
	    int lt14 = lminpo;
	    while (lt14 <= lmaxpo) {
	      i3j0++;
	      int l3 = lt14 - 1;
	      int ny = l3 * l3 + lt14;
	      double aors = 1.0 * (l3 + lt14);
	      double f3j = (c1->v3j0[i3j0 - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	      cfun = (c1->vh[nbhj + lt14 - 1] * c1->vyhj[nby + ny - 1]) * f3j;
	      csum += cfun;
	      jsn *= -1;
	      lt14 += 2;
	    }
	    // goes to 22
	  } else { // label 16
	    r3jjr(l1mp, l2, -mupm1, mupm2, c1->rac3j);
	    int il = ilin;
	    int lt20 = lminpo;
	    while (lt20 <= lmaxpo) {
	      i3j0++;
	      if (m1mm2m <= lt20) {
		il += 2;
		int l3 = lt20 - 1;
		int ny = l3 * l3  + lt20 + m1mm2;
		double aors = 1.0 * (l3 + lt20);
		double f3j = (c1->rac3j[il - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
		cfun = (c1->vh[nbhj + lt20 - 1] * c1->vyhj[nby + ny - 1]) * f3j;
		csum += cfun;
	      }
	      // label 20
	      jsn *= -1;
	      lt20 += 2;
	    }
	  }
	  // label 22
	  csum *= cr;
	  result += csum;
	}
	// Otherwise there is nothing to add
      } // jm24 loop. Should go to 70
    } else { // label 30, IHI == 1
      for (int jm44 = 1; jm44 <= 3; jm44++) {
	csum = cc0;
	int mu = jm44 - 2;
	int mupm1 = mu + m1;
	int mupm2 = mu + m2;
	if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	  int jsn = - isn;
	  if (mu == 0) jsn = isn;
	  double cr = cgev(ipamo, mu, l1, m1) * cgev(0, mu, l2, m2);
	  int i3j0 = i3j0in;
	  if (mupm1 == 0 && mupm2 == 0) {
	    int lt34 = lminpo;
	    while (lt34 <= lmaxpo) {
	      i3j0++;
	      int l3 = lt34 - 1;
	      int ny = l3 * l3 + lt34;
	      double aors = 1.0 * (l3 + lt34);
	      double f3j = (c1->v3j0[i3j0 - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	      cfun = (c1->vj * c1->vyhj[nby + ny - 1]) * f3j;
	      csum += cfun;
	      jsn *= -1;
	      lt34 += 2;
	    }
	    // goes to 42
	  } else { // label 36
	    r3jjr(l1mp, l2, -mupm1, mupm2, c1->rac3j);
	    int il = ilin;
	    int lt40 = lminpo;
	    while (lt40 <= lmaxpo) {
	      i3j0++;
	      if (m1mm2m <= lt40) {
		il += 2;
		int l3 = lt40 - 1;
		int ny = l3 * l3  + lt40 + m1mm2;
		double aors = 1.0 * (l3 + lt40);
		double f3j = (c1->rac3j[il - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
		cfun = (c1->vj * c1->vyhj[nby + ny - 1]) * f3j;
		csum += cfun;
	      }
	      // label 40
	      jsn *= -1;
	      lt40 += 2;
	    }
	  }
	  // label 42
	  csum *= cr;
	  result += csum;
	}
	// Otherwise there is nothing to add
      } // jm44 loop. Should go to 70
    }
    // goes to 70
  } else { // label 50, IHI == 2
    int nbhj = nblmo * c1->lmtpo;
    int nby = nblmo * c1->lmtpos;
    for (int jm64 = 1; jm64 <= 3; jm64++) {
      csum = cc0;
      int mu = jm64 - 2;
      int mupm1 = mu + m1;
      int mupm2 = mu + m2;
      if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	int jsn = -isn;
	if (mu == 0) jsn = isn;
	double cr = cgev(ipamo, mu, l1, m1) * cgev(0, mu, l2, m2);
	int i3j0 = i3j0in;
	if (mupm1 == 0 && mupm2 == 0) {
	  int lt54 = lminpo;
	  while (lt54 <= lmaxpo) {
	    i3j0++;
	    int l3 = lt54 - 1;
	    int ny = l3 * l3 + lt54;
	    double aors = 1.0 * (l3 + lt54);
	    double f3j = (c1->v3j0[i3j0 - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	    cfun = (c1->vj0[nbhj + lt54 - 1] * c1->vyj0[nby + ny - 1]) * f3j;
	    csum += cfun;
	    jsn *= -1;
	    lt54 += 2;
	  }
	  // goes to 62
	} else { // label 56
	  r3jjr(l1mp, l2, -mupm1, mupm2, c1->rac3j);
	  int il = ilin;
	  int lt60 = lminpo;
	  while (lt60 <= lmaxpo) {
	    i3j0++;
	    if (m1mm2m <= lt60) {
	      il += 2;
	      int l3 = lt60 - 1;
	      int ny = l3 * l3  + lt60 + m1mm2;
	      double aors = 1.0 * (l3 + lt60);
	      double f3j = (c1->rac3j[il - 1] * c1->v3j0[i3j0 - 1] * sqrt(aors)) * jsn;
	      cfun = (c1->vj0[nbhj + lt60 - 1] * c1->vyj0[nby + ny - 1]) * f3j;
	      csum += cfun;
	    }
	    // label 60
	    jsn *= -1;
	    lt60 += 2;
	  }
	}
	// label 62
	csum *= cr;
	result += csum;
      }
      // Otherwise there is nothing to add
    } // jm64 loop. Should go to 70
  }
  // label 70
  const double four_pi = acos(0.0) * 8.0;
  if (ipamo != 1) {
    double cr = sqrt(four_pi * (l1 + l1po) * (l2 + l2 + 1));
    result *= cr;
  } else {
    double cr = sqrt(four_pi * (l1 + l1mp) * (l1 + l1po) * (l2 + l2 + 1) / l1po);
    result *= (cr * uim);
  }
  return result;
}
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp end declare target
// #endif

void hjv(
	 double exri, double vk, int &jer, int &lcalc, dcomplex &arg,
	 ParticleDescriptor *c1
) {
  int nsphmo = c1->nsph - 1;
  int lit = c1->li + c1->li;
  int lmt = c1->li + c1->le;
  const int rfj_size = (lit > lmt) ? lit : lmt;
  const int rfn_size = c1->litpo;
  double *rfj, *rfn;
  rfj = new double[rfj_size+1]();
  rfn = new double[rfn_size+1]();
  jer = 0;
  int ivhb = 0;
  for (int nf40 = 1; nf40 <= nsphmo; nf40++) { // GPU portable?
    int nfpo = nf40 + 1;
    for (int ns40 = nfpo; ns40 <= c1->nsph; ns40++) {
      double rx = c1->rxx[nf40 - 1] - c1->rxx[ns40 - 1];
      double ry = c1->ryy[nf40 - 1] - c1->ryy[ns40 - 1];
      double rz = c1->rzz[nf40 - 1] - c1->rzz[ns40 - 1];
      double rr = sqrt(rx * rx + ry * ry + rz * rz);
      double rarg = rr * vk * exri;
      arg = rarg +  0.0 * I;
      rbf(lit, rarg, lcalc, rfj);
      if (lcalc < lit) {
	jer = 1;
	delete[] rfj;
	delete[] rfn;
	return;
      }
      rnf(lit, rarg, lcalc, rfn);
      if (lcalc < lit) {
	jer = 2;
	delete[] rfj;
	delete[] rfn;
	return;
      }
      for (int lpo38 = 1; lpo38 <= c1->litpo; lpo38++) {
	double rpart = rfj[lpo38 - 1];
	double ipart = rfn[lpo38 - 1];
	c1->vh[lpo38 + ivhb - 1] = rpart + ipart * I;
      }
      ivhb += c1->litpo;
    } // ns40 loop
  } // nf40 loop
  ivhb = 0;
  for (int nf50 = 1; nf50 <= c1->nsph; nf50++) {
    double rx = c1->rxx[nf50 - 1];
    double ry = c1->ryy[nf50 - 1];
    double rz = c1->rzz[nf50 - 1];
    if (!(rx == 0.0 && ry == 0.0 && rz == 0.0)) {
      double rr = sqrt(rx * rx + ry * ry + rz * rz);
      double rarg = rr * vk * exri;
      rbf(lmt, rarg, lcalc, rfj);
      if (lcalc < lmt) {
	jer = 3;
	delete[] rfj;
	delete[] rfn;
	return;
      }
      for (int lpo47 = 1; lpo47 <= c1->lmtpo; lpo47++) {
	c1->vj0[lpo47 + ivhb - 1] = rfj[lpo47 - 1];
      }
    }
    ivhb += c1->lmtpo;
  } // nf50 loop
  delete[] rfj;
  delete[] rfn;
}

void lucin(dcomplex **am, const np_int nddmst, np_int n, int &ier) {
  /* NDDMST  FIRST DIMENSION OF AM AS DECLARED IN DIMENSION
   *         STATEMENT.
   * N       NUMBER OF ROWS IN AM.
   * IER     IS REPLACED BY 1 FOR SINGULARITY.
   */
  double *v = new double[nddmst];
  dcomplex ctemp, cfun;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  ier = 0;
  int nminus = n - 1;
  for (int64_t i = 1; i <= n; i++) {
    double sum = 0.0;
    for (int64_t j = 1; j <= n; j++) {
      sum += (
	      real(am[i - 1][j - 1]) * real(am[i - 1][j - 1])
	      + imag(am[i - 1][j - 1]) * imag(am[i - 1][j - 1])
	      );
    } // j1319 loop
    v[i - 1] = 1.0 / sum;
  } // i1309 loop
  // 2. REPLACE AM BY TRIANGULAR MATRICES (L,U) WHERE AM=L*U.
  //    REPLACE L(I,I) BY 1/L(I,I), READY FOR SECTION 4.
  //    (ROW INTERCHANGES TAKE PLACE, AND THE INDICES OF THE PIVOTAL ROWS
  //    ARE PLACED IN V.)
  for (int64_t k = 1; k <= n; k++) {
    int64_t kplus = k + 1;
    int64_t kminus = k - 1;
    int64_t l = k;
    double psqmax = 0.0;
    for (int64_t i = k; i <= n; i++) {
      cfun = cdtp(-am[i - 1][k - 1], am, i, 1, k, kminus);
      ctemp = -cfun;
      am[i - 1][k - 1] = ctemp;
      double psq = v[i - 1] * (real(ctemp) * real(ctemp) + imag(ctemp) * imag(ctemp));
      if (psq > psqmax) {
	psqmax = psq;
	l = i;
      }
    } // i2029 loop
    if (l != k) {
      for (int64_t j = 1; j <= n; j++) {
	ctemp = am[k - 1][j - 1];
	am[k - 1][j - 1] = am[l - 1][j - 1];
	am[l - 1][j - 1] = ctemp;
      } // j2049 loop
      v[l - 1] = v[k - 1];
    }
    // label 2011
    v[k - 1] = 1.0 * l;
    if (psqmax == 0.0) {
      ier = 1;
      delete[] v;
      return;
    }
    ctemp = 1.0 / am[k - 1][k - 1];
    am[k - 1][k - 1] = ctemp;
    if (kplus <= n) {
      for (int64_t j = kplus; j <= n; j++) {
	cfun = cdtp(-am[k - 1][j - 1], am, k, 1, j, kminus);
	am[k - 1][j - 1] = -ctemp * cfun;
      } // j2059 loop
    }
  } // k2019 loop
  // 4.  REPLACE AM BY ITS INVERSE AMINV.
  // 4.1 REPLACE L AND U BY THEIR INVERSE LINV AND UINV.
  for (int64_t k = 1; k <= nminus; k++) {
    int64_t kplus = k + 1;
    for (int64_t i = kplus; i <= n; i++) {
      cfun = cdtp(cc0, am, i, k, k, i - k);
      am[i - 1][k - 1] = -am[i - 1][i - 1] * cfun;
      cfun = cdtp(am[k - 1][i - 1], am, k, kplus, i, i - k - 1);
      am[k - 1][i - 1] = -cfun;
    } // i4119 loop
  } // k4109 loop
  // 4.2 FORM AMINV=UINV*LINV.
  for (int64_t k = 1; k <= n; k++) {
    for (int64_t i = 1; i <= n; i++) {
      if (i < k) {
	cfun = cdtp(cc0, am, i, k, k, n - k + 1);
	am[i - 1][k -1] = cfun;
      }
      else {
	cfun = cdtp(am[i - 1][k - 1], am, i, i + 1, k, n - i);
	am[i - 1][k - 1] = cfun;
      }
    } // i4119 loop
  } // k4209 loop
  // 4.3 INTERCHANGE COLUMNS OF AMINV AS SPECIFIED BY V, BUT IN REVERSE
  //     ORDER.
  for (int64_t l = 1; l <= n; l++) {
    int64_t k = n - l + 1;
    int64_t kcol = (int64_t)(v[k - 1]);
    if (kcol != k) {
      for (int64_t i = 1; i <= n; i++) {
	ctemp = am[i - 1][k - 1];
	am[i - 1][k - 1] = am[i - 1][kcol - 1];
	am[i - 1][kcol - 1] = ctemp;
      } // i4319 loop
    }
  } // l4309 loop
  delete[] v;
}

void mextc(double vk, double exri, dcomplex **fsac, double **cextlr, double **cext) {
  double fa11r = real(fsac[0][0]);
  double fa11i = imag(fsac[0][0]);
  double fa21r = real(fsac[1][0]);
  double fa21i = imag(fsac[1][0]);
  double fa12r = real(fsac[0][1]);
  double fa12i = imag(fsac[0][1]);
  double fa22r = real(fsac[1][1]);
  double fa22i = imag(fsac[1][1]);
  cextlr[0][0] = fa11i * 2.0;
  cextlr[0][1] = 0.0;
  cextlr[0][2] = -fa12i;
  cextlr[0][3] = -fa12r;
  cextlr[1][0] = 0.0;
  cextlr[1][1] = fa22i * 2.0;
  cextlr[1][2] = -fa21i;
  cextlr[1][3] = fa21r;
  cextlr[2][0] = -fa21i * 2.0;
  cextlr[2][1] = -fa12i * 2.0;
  cextlr[2][2] = fa11i + fa22i;
  cextlr[2][3] = fa22r - fa11r;
  cextlr[3][0] = fa21r * 2.0;
  cextlr[3][1] = -fa12r * 2.0;
  cextlr[3][2] = fa11r - fa22r;
  cextlr[3][3] = cextlr[2][2];
  cext[0][0] = cextlr[3][3];
  cext[1][1] = cextlr[3][3];
  cext[2][2] = cextlr[3][3];
  cext[2][3] = cextlr[2][3];
  cext[3][2] = cextlr[3][2];
  cext[3][3] = cextlr[3][3];
  cext[0][1] = fa11i - fa22i;
  cext[0][2] = -fa12i - fa21i;
  cext[0][3] = fa21r - fa12r;
  cext[1][0] = cext[0][1];
  cext[1][2] = fa21i - fa12i;
  cext[3][1] = fa12r + fa21r;
  cext[1][3] = -cext[3][1];
  cext[2][0] = cext[0][2];
  cext[2][1] = -cext[1][2];
  cext[3][0] = cext[1][3];
  double ckm = vk / exri;
  for (int i10 = 0; i10 < 4; i10++) {
    for (int j10 = 0; j10 < 4; j10++) {
      cextlr[i10][j10] *= ckm;
      cext[i10][j10] *= ckm;
    }
  }
}

void pcros(double vk, double exri, ParticleDescriptor *c1) {
#ifdef USE_NVTX
  nvtxRangePush("whole pcros");
#endif
  const dcomplex cc0 = 0.0 + 0.0 * I;
  dcomplex sump, sum1, sum2, sum3, sum4, am, amp, cc, csam;
  const double exdc = exri * exri;
  double ccs = 1.0 / (vk * vk);
  double cccs = ccs / exdc;
  csam = -(ccs / (exri * vk)) * 0.5 * I;
  const double pi4sq = 64.0 * acos(0.0) * acos(0.0);
  double cfsq = 4.0 / (pi4sq *ccs * ccs);
  const int nlemt = c1->nlem + c1->nlem;
  int jpo = 2;
  dcomplex *vec_am0m = c1->am0m[0];
  dcomplex *vec_w = c1->w[0];
#ifdef USE_NVTX
  nvtxRangePush("pcros outer loop 1");
#endif
  for (int ipo18 = 0; ipo18 < 2; ipo18++) {
    int jpo = 1-ipo18;
    int ipopt = ipo18 + 2;
    int jpopt = jpo + 2;
    double sum = 0.0;
    sump = cc0;
    sum1 = cc0;
    sum2 = cc0;
    sum3 = cc0;
    sum4 = cc0;
#ifdef USE_NVTX
  nvtxRangePush("pcros intermediate loop 1");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:sum, sump, sum1, sum2, sum3, sum4)
// #else
// #pragma omp parallel for simd reduction(+:sum, sump, sum1, sum2, sum3, sum4)
// #endif
#pragma omp parallel for simd reduction(+:sum, sump, sum1, sum2, sum3, sum4)
  for (int i12 = 0; i12 < nlemt; i12++) {
      // int i = i12 - 1;
      dcomplex am = cc0;
      dcomplex amp = cc0;
      for (int j10 = 0; j10 < nlemt; j10++) {
	// int j = j10 - 1;
	am += (vec_am0m[nlemt*i12+j10] * vec_w[4*j10+ipo18]);
	amp += (vec_am0m[nlemt*i12+j10] * vec_w[4*j10+jpo]);
      } // j10 loop
      sum += real(dconjg(am) * am);
      sump += (dconjg(amp) * am);
      sum1 += (dconjg(vec_w[4*i12+ipo18]) * am);
      sum2 += (dconjg(vec_w[4*i12+jpo]) * am);
      sum3 += (vec_w[4*i12+ipopt] * am);
      sum4 += (vec_w[4*i12+jpopt] * am);
    } // i12 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
    c1->scsc[ipo18] = cccs * sum;
    c1->scscp[ipo18] = cccs * sump;
    c1->ecsc[ipo18] = -cccs * real(sum1);
    c1->ecscp[ipo18] = -cccs * sum2;
    c1->fsac[ipo18][ipo18] = csam * sum1;
    c1->fsac[jpo][ipo18] = csam * sum2;
    c1->sac[ipo18][ipo18] = csam * sum3;
    c1->sac[jpo][ipo18] = csam * sum4;
  } // ipo18 loop
  int i = 0;
  dcomplex * &vint =  c1->vint;
#ifdef USE_NVTX
  nvtxRangePush("pcros loop 2");
#endif
  for (int ipo1 = 1; ipo1 <= 2; ipo1++) {
    for (int jpo1 = 1; jpo1 <= 2; jpo1++) {
      cc = dconjg(c1->sac[jpo1 - 1][ipo1 - 1]);
      for (int ipo2 = 1; ipo2 <= 2; ipo2 ++) {
	for (int jpo2 = 1; jpo2 <= 2; jpo2++) {
	  c1->vint[i++] = c1->sac[jpo2 - 1][ipo2 - 1] * cc * cfsq;
	} // jpo2 loop
      } // ipo2 loop
    } // jpo1 loop
  } // ipo1 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}

void pcrsm0(double vk, double exri, int inpol, ParticleDescriptor *c1) {
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex sum1, sum2, sum3, sum4;
  dcomplex sums1, sums2, sums3, sums4, csam;
  double exdc = exri * exri;
  double ccs = 4.0 * acos(0.0) / (vk * vk);
  double cccs = ccs / exdc;
  int nlemt = c1->nlem + c1->nlem;
  dcomplex *vec_am0m = c1->am0m[0];
  csam = -(ccs / (exri * vk)) * 0.5 * I;
  sum2 = cc0;
  sum3 = cc0;
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:sum2,sum3)
// #else
// #pragma omp parallel for simd reduction(+:sum2,sum3)
// #endif
#pragma omp parallel for simd reduction(+:sum2,sum3)
  for (int i14 = 0; i14 < c1->nlem; i14++) { 
    int ie = i14 + c1->nlem;
    sum2 += (vec_am0m[nlemt*i14 + i14] + vec_am0m[nlemt*ie + ie]);
    sum3 += (vec_am0m[nlemt*i14 + ie] + vec_am0m[nlemt*ie + i14]);
  } // i14 loop
  double sumpi = 0.0;
  dcomplex sumpd = cc0;
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd collapse(2) reduction(+:sumpi,sumpd)
// #else
// #pragma omp parallel for simd collapse(2) reduction(+:sumpi,sumpd)
// #endif
#pragma omp parallel for simd collapse(2) reduction(+:sumpi,sumpd)
  for (int i16 = 0; i16 < nlemt; i16++) {
    for (int j16 = 0; j16 < c1->nlem; j16++) {
      int je = j16 + c1->nlem;
      double rvalue = real(
			   dconjg(vec_am0m[nlemt*i16 + j16]) * vec_am0m[nlemt*i16 + j16]
			   + dconjg(vec_am0m[nlemt*i16 + je]) * vec_am0m[nlemt*i16 + je]
			   );
      sumpi += rvalue;
      sumpd += (
		dconjg(vec_am0m[nlemt*i16 + j16]) * vec_am0m[nlemt*i16 + je]
		+ dconjg(vec_am0m[nlemt*i16 + je]) * vec_am0m[nlemt*i16 + j16]
		);
    } // j16 loop
  } // i16 loop
  if (inpol == 0) {
    sum1 = sum2;
    sum4 = sum3 * uim;
    sum3 = -sum4;
    sums1 = sumpi;
    sums2 = sumpi;
    sums3 = sumpd * uim;
    sums4 = -sums3;
  } else { // label 18
    sum1 = sum2 + sum3;
    sum2 = sum2 - sum3;
    sum3 = cc0;
    sum4 = cc0;
    sums1 = sumpi - sumpd;
    sums2 = sumpi + sumpd;
    sums3 = cc0;
    sums4 = cc0;
  }
  // label 20
  c1->ecscm[0] = -cccs * real(sum2);
  c1->ecscm[1] = -cccs * real(sum1);
  c1->ecscpm[0] = -cccs * sum4;
  c1->ecscpm[1] = -cccs * sum3;
  c1->fsacm[0][0] = csam * sum2;
  c1->fsacm[1][0] = csam * sum4;
  c1->fsacm[1][1] = csam * sum1;
  c1->fsacm[0][1] = csam * sum3;
  c1->scscm[0] = cccs * real(sums1);
  c1->scscm[1] = cccs * real(sums2);
  c1->scscpm[0] = cccs * sums3;
  c1->scscpm[1] = cccs * sums4;
}

void polar(
	   double x, double y, double z, double &r, double &cth, double &sth,
	   double &cph, double &sph
) {
  bool onx = (y == 0.0);
  bool ony = (x == 0.0);
  bool onz = (onx && ony);
  double rhos = 0.0;
  double rho = 0.0;
  if (!onz) {
    if (!onx) {
      if (!ony) {
	rhos = x * x + y * y;
	rho = sqrt(rhos);
	cph = x / rho;
	sph = y / rho;
	// goes to 25
      } else { // label 20
	rhos = y * y;
	rho = (y > 0.0) ? y : -y;
	cph = 0.0;
	sph = (y > 0.0) ? 1.0 : -1.0;
	// goes to 25
      }
    } else { // label 15
      rhos = x * x;
      rho = (x > 0.0) ? x : -x;
      cph = (x > 0.0) ? 1.0 : -1.0;
      sph = 0.0;
      // goes to 25
    }
  } else { // label 10
    cph = 1.0;
    sph = 0.0;
    // goes to 25
  }
  // label 25
  if (z == 0.0) {
    if (!onz) {
      r = rho;
      cth = 0.0;
      sth = 1.0;
      // returns
    } else { // label 30
      r = 0.0;
      cth = 1.0;
      sth = 0.0;
      // returns
    }
  } else { // label 35
    if (!onz) {
      r = sqrt(rhos + z * z);
      cth = z / r;
      sth = rho / r;
      // returns
    } else { // label 40
      r = (z > 0.0) ? z : -z;
      cth = (z > 0.0) ? 1.0 : -1.0;
      sth = 0.0;
      // returns
    }
  }
}

void r3j000(int j2, int j3, double *rac3j) {
  int jmx = j3 + j2;
  if (jmx <= 0) {
    rac3j[0] = 1.0;
    return;
  }
  int jmn = j3 - j2;
  if (jmn < 0) jmn *= -1;
  int njmo = (jmx - jmn) / 2;
  int jf = jmx + jmx + 1;
  int isn = 1;
  if (jmn % 2 != 0) isn = -1;
  if (njmo <= 0) {
    double sj = 1.0 * jf;
    double cnr = (1 / sqrt(sj)) * isn;
    rac3j[0] = cnr;
    return;
  }
  double sjr = 1.0 * jf;
  int jmxpos = (jmx + 1) * (jmx + 1);
  int jmns = jmn * jmn;
  int j1mo = jmx - 1;
  int j1s = (j1mo + 1) * (j1mo + 1);
  double cj = sqrt(1.0 * (jmxpos - j1s) * (j1s - jmns));
  int j1mos = j1mo * j1mo;
  double cjmo = sqrt(1.0 * (jmxpos - j1mos) * (j1mos - jmns));
  if (njmo <= 1) {
    rac3j[0] = -cj / cjmo;
    double sj = sjr + (rac3j[0] * rac3j[0]) * (jf - 4);
    double cnr = (1.0 / sqrt(sj)) * isn;
    rac3j[1] = cnr;
    rac3j[0] *= cnr;
    return;
  }
  int nj = njmo + 1;
  int nmat = (nj + 1) / 2;
  rac3j[nj - 1] = 1.0;
  rac3j[njmo - 1] = -cj / cjmo;
  if (nmat != njmo) {
    int nbr = njmo - nmat;
    for (int ibr45 = 1; ibr45 <= nbr; ibr45++) {
      int irr = nj - ibr45;
      jf -= 4;
      j1mo -= 2;
      j1s = (j1mo + 1) * (j1mo + 1);
      cj = sqrt(1.0 * (jmxpos - j1s) * (j1s - jmns));
      j1mos = j1mo * j1mo;
      cjmo = sqrt(1.0 * (jmxpos - j1mos) * (j1mos - jmns));
      rac3j[irr - 2] = rac3j[irr - 1] * (-cj / cjmo);
      sjr = sjr + (rac3j[irr - 1] * rac3j[irr - 1]) * jf;
    }
  }
  // label 50
  double racmat = rac3j[nmat - 1];
  sjr = sjr + (racmat * racmat) * (jf - 4);
  rac3j[0] = 1.0;
  jf = jmn + jmn + 1;
  double sjl = 1.0 * jf;
  int j1pt = jmn + 2;
  int j1pos = (j1pt - 1) * (j1pt - 1);
  double cjpo = sqrt(1.0 * (jmxpos - j1pos) * (j1pos - jmns));
  int j1pts = j1pt * j1pt;
  double cjpt = sqrt(1.0 * (jmxpos - j1pts) * (j1pts - jmns));
  rac3j[1] = -cjpo / cjpt;
  int nmatmo = nmat - 1;
  if (nmatmo >= 2) {
    for (int irl70 = 2; irl70 <= nmatmo; irl70++) {
      jf += 4;
      j1pt += 2;
      j1pos = (j1pt - 1) * (j1pt - 1);
      cjpo = sqrt(1.0 * (jmxpos - j1pos) * (j1pos - jmns));
      j1pts = j1pt * j1pt;
      cjpt = sqrt(1.0 * (jmxpos - j1pts) * (j1pts - jmns));
      rac3j[irl70] = rac3j[irl70 - 1] * (-cjpo / cjpt);
      sjl = sjl + (rac3j[irl70 - 1] * rac3j[irl70 - 1]) * jf;
    }
  }
  // label 75
  double ratrac = racmat / rac3j[nmat - 1];
  double rats = ratrac * ratrac;
  double sj = sjr + sjl * rats;
  rac3j[nmat - 1] = racmat;
  double cnr = (1.0 / sqrt(sj)) * isn;
  for (int irr80 = nmat; irr80 <= nj; irr80++) {
    rac3j[irr80 - 1] *= cnr;
  }
  double cnl = cnr * ratrac;
  for (int irl85 = 1; irl85 <= nmatmo; irl85++) {
    rac3j[irl85 - 1] *= cnl;
  }
}

// #ifdef USE_TARGET_OFFLOAD
// #pragma omp begin declare target device_type(any)
// #endif
void r3jjr(int j2, int j3, int m2, int m3, double *rac3j) {
  int jmx = j3 + j2;
  int jdf = j3 - j2;
  int m1 = -m2 - m3;
  int abs_jdf = (jdf >= 0) ? jdf : -jdf;
  int abs_m1 = (m1 >= 0) ? m1 : -m1;
  int jmn = (abs_jdf > abs_m1) ? abs_jdf : abs_m1;
  int njmo = jmx - jmn;
  int jf = jmx + jmx + 1;
  int isn = 1;
  if ((jdf + m1) % 2 != 0) isn = -1;
  if (njmo <= 0) {
    double sj = 1.0 * jf;
    double cnr = (1.0 / sqrt(sj)) * isn;
    rac3j[0] = cnr;
  } else { // label 15
    double sjt = 1.0;
    double sjr = 1.0 * jf;
    int jsmpos = (jmx + 1) * (jmx + 1);
    int jdfs = jdf * jdf;
    int m1s = m1 * m1;
    int mdf = m3 - m2;
    int idjc = m1 * (j3 * (j3 + 1) - j2 * (j2 +1));
    int j1 = jmx;
    int j1s = j1 * j1;
    int j1po = j1 + 1;
    double ccj = 1.0 * (j1s - jdfs) * (j1s - m1s);
    double cj = sqrt(ccj * (jsmpos - j1s));
    double dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
    if (njmo <= 1) {
      rac3j[0] = -dj / (cj * j1po);
      double sj = sjr + (rac3j[0] * rac3j[0]) * (jf - 2);
      double cnr = (1.0 / sqrt(sj)) * isn;
      rac3j[1] = cnr;
      rac3j[0] *= cnr;
    } else { // label 20
      double cjp = 0.0;
      int nj = njmo + 1;
      int nmat = (nj + 1) / 2;
      rac3j[nj - 1] = 1.0;
      rac3j[njmo - 1] = -dj / (cj * j1po);
      if (nmat != njmo) {
	int nbr = njmo - nmat;
	for (int ibr45 = 1; ibr45 <= nbr; ibr45++) {
	  int irr = nj - ibr45;
	  jf -= 2;
	  j1--;
	  j1s = j1 * j1;
	  j1po = j1 + 1;
	  cjp = cj;
	  ccj = 1.0 * (j1s - jdfs) * (j1s - m1s);
	  cj = sqrt(ccj * (jsmpos - j1s));
	  sjt = rac3j[irr - 1] * rac3j[irr - 1];
	  dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
	  rac3j[irr - 2] = -(rac3j[irr - 1] * dj
				 + rac3j[irr] * cjp * j1) / (cj * j1po);
	  sjr += (sjt * jf);
	} // ibr45 loop
      }
      // label 50
      double osjt = sjt;
      sjt = rac3j[nmat - 1] * rac3j[nmat - 1];
      if (sjt >= osjt) {
	sjr += (sjt * (jf - 2));
      } else { // label 55
	nmat++;
      }
      // label 60
      double racmat = rac3j[nmat - 1];
      rac3j[0] = 1.0;
      jf = jmn + jmn + 1;
      double sjl = 1.0 * jf;
      j1 = jmn;
      if (j1 != 0) {
	j1po = j1 + 1;
	int j1pos = j1po * j1po;
	double ccjp = 1.0 * (j1pos - jdfs) * (j1pos - m1s);
	cjp = sqrt(ccjp * (jsmpos - j1pos));
	dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
	rac3j[1] = - dj / (cjp * j1);
      } else { // label 62
	cjp = sqrt(1.0 * (jsmpos - 1));
	dj = 1.0 * mdf;
	rac3j[1] = -dj / cjp;
      }
      // label 63
      int nmatmo = nmat - 1;
      if (nmatmo >= 2) {
	for (int irl70 = 2; irl70 <= nmatmo; irl70++) {
	  jf += 2;
	  j1++;
	  j1po = j1 + 1;
	  int j1pos = j1po * j1po;
	  cj = cjp;
	  double ccjp = 1.0 * (j1pos - jdfs) * (j1pos - m1s);
	  cjp = sqrt(ccjp * (jsmpos - j1pos));
	  sjt = rac3j[irl70 - 1] * rac3j[irl70 - 1];
	  dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
	  rac3j[irl70] = -(
			   rac3j[irl70 - 1] * dj
			   + rac3j[irl70 - 2] * cj * j1po
			   ) / (cjp * j1);
	  sjl += (sjt * jf);
	}
      }
      // label 75
      double ratrac = racmat / rac3j[nmat - 1];
      double rats = ratrac * ratrac;
      double sj = sjr + sjl * rats;
      rac3j[nmat - 1] = racmat;
      double cnr = (1.0 / sqrt(sj)) * isn;
      for (int irr80 = nmat; irr80 <= nj; irr80++) rac3j[irr80 - 1] *= cnr;
      double cnl = cnr * ratrac;
      for (int irl85 = 1; irl85 <= nmatmo; irl85++) rac3j[irl85 - 1] *= cnl;
    }
  }
}
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp end declare target
// #endif

// #ifdef USE_TARGET_OFFLOAD
// #pragma omp begin declare target device_type(any)
// #endif
void r3jjr_d(int j2, int j3, int m2, int m3, double *rac3j) {
  int jmx = j3 + j2;
  int jdf = j3 - j2;
  int m1 = -m2 - m3;
  int abs_jdf = (jdf >= 0) ? jdf : -jdf;
  int abs_m1 = (m1 >= 0) ? m1 : -m1;
  int jmn = (abs_jdf > abs_m1) ? abs_jdf : abs_m1;
  int njmo = jmx - jmn;
  int jf = jmx + jmx + 1;
  int isn = 1;
  if ((jdf + m1) % 2 != 0) isn = -1;
  if (njmo <= 0) {
    double sj = 1.0 * jf;
    double cnr = (1.0 / sqrt(sj)) * isn;
    rac3j[0] = cnr;
  } else { // label 15
    double sjt = 1.0;
    double sjr = 1.0 * jf;
    int jsmpos = (jmx + 1) * (jmx + 1);
    int jdfs = jdf * jdf;
    int m1s = m1 * m1;
    int mdf = m3 - m2;
    int idjc = m1 * (j3 * (j3 + 1) - j2 * (j2 +1));
    int j1 = jmx;
    int j1s = j1 * j1;
    int j1po = j1 + 1;
    double ccj = 1.0 * (j1s - jdfs) * (j1s - m1s);
    double cj = sqrt(ccj * (jsmpos - j1s));
    double dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
    if (njmo <= 1) {
      rac3j[0] = -dj / (cj * j1po);
      double sj = sjr + (rac3j[0] * rac3j[0]) * (jf - 2);
      double cnr = (1.0 / sqrt(sj)) * isn;
      rac3j[1] = cnr;
      rac3j[0] *= cnr;
    } else { // label 20
      double cjp = 0.0;
      int nj = njmo + 1;
      int nmat = (nj + 1) / 2;
      rac3j[nj - 1] = 1.0;
      rac3j[njmo - 1] = -dj / (cj * j1po);
      if (nmat != njmo) {
	int nbr = njmo - nmat;
	for (int ibr45 = 1; ibr45 <= nbr; ibr45++) {
	  int irr = nj - ibr45;
	  jf -= 2;
	  j1--;
	  j1s = j1 * j1;
	  j1po = j1 + 1;
	  cjp = cj;
	  ccj = 1.0 * (j1s - jdfs) * (j1s - m1s);
	  cj = sqrt(ccj * (jsmpos - j1s));
	  sjt = rac3j[irr - 1] * rac3j[irr - 1];
	  dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
	  rac3j[irr - 2] = -(rac3j[irr - 1] * dj
				 + rac3j[irr] * cjp * j1) / (cj * j1po);
	  sjr += (sjt * jf);
	} // ibr45 loop
      }
      // label 50
      double osjt = sjt;
      sjt = rac3j[nmat - 1] * rac3j[nmat - 1];
      if (sjt >= osjt) {
	sjr += (sjt * (jf - 2));
      } else { // label 55
	nmat++;
      }
      // label 60
      double racmat = rac3j[nmat - 1];
      rac3j[0] = 1.0;
      jf = jmn + jmn + 1;
      double sjl = 1.0 * jf;
      j1 = jmn;
      if (j1 != 0) {
	j1po = j1 + 1;
	int j1pos = j1po * j1po;
	double ccjp = 1.0 * (j1pos - jdfs) * (j1pos - m1s);
	cjp = sqrt(ccjp * (jsmpos - j1pos));
	dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
	rac3j[1] = - dj / (cjp * j1);
      } else { // label 62
	cjp = sqrt(1.0 * (jsmpos - 1));
	dj = 1.0 * mdf;
	rac3j[1] = -dj / cjp;
      }
      // label 63
      int nmatmo = nmat - 1;
      if (nmatmo >= 2) {
	for (int irl70 = 2; irl70 <= nmatmo; irl70++) {
	  jf += 2;
	  j1++;
	  j1po = j1 + 1;
	  int j1pos = j1po * j1po;
	  cj = cjp;
	  double ccjp = 1.0 * (j1pos - jdfs) * (j1pos - m1s);
	  cjp = sqrt(ccjp * (jsmpos - j1pos));
	  sjt = rac3j[irl70 - 1] * rac3j[irl70 - 1];
	  dj = 1.0 * jf * (j1 * j1po * mdf + idjc);
	  rac3j[irl70] = -(
			       rac3j[irl70 - 1] * dj
			       + rac3j[irl70 - 2] * cj * j1po
			       ) / (cjp * j1);
	  sjl += (sjt * jf);
	}
      }
      // label 75
      double ratrac = racmat / rac3j[nmat - 1];
      double rats = ratrac * ratrac;
      double sj = sjr + sjl * rats;
      rac3j[nmat - 1] = racmat;
      double cnr = (1.0 / sqrt(sj)) * isn;
      for (int irr80 = nmat; irr80 <= nj; irr80++) rac3j[irr80 - 1] *= cnr;
      double cnl = cnr * ratrac;
      for (int irl85 = 1; irl85 <= nmatmo; irl85++) rac3j[irl85 - 1] *= cnl;
    }
  }
}
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp end declare target
// #endif

void r3jmr(int j1, int j2, int j3, int m1, double *rac3j) {
  int mmx = (j2 < j3 - m1) ? j2 : j3 - m1;
  int mmn = (-j2 > -(j3 + m1)) ? -j2 : -(j3 + m1);
  int nmmo = mmx - mmn;
  int j1po = j1 + 1;
  int j1tpo = j1po + j1;
  int isn = 1;
  if ((j2 - j3 - m1) % 2 != 0) isn = -1;
  if (nmmo <= 0) {
    double sj = 1.0 * j1tpo;
    double cnr = (1.0 / sqrt(sj)) * isn;
    rac3j[0] = cnr;
    // returns
  } else { // label 15
    int j1s = j1 * j1po;
    int j2po = j2 + 1;
    int j2s = j2 * j2po;
    int j3po = j3 + 1;
    int j3s = j3 * j3po;
    int id = j1s - j2s - j3s;
    int m2 = mmx;
    int m3 = m1 + m2;
    double cm = sqrt(1.0 * (j2po - m2) * (j2 + m2) * (j3po - m3) * (j3 + m3));
    double dm = 1.0 * (id + m2 * m3 * 2);
    if (nmmo <= 1) {
      rac3j[0] = dm / cm;
      double sj = (1.0 + rac3j[0] * rac3j[0]) * j1tpo;
      double cnr = 1.0 / sqrt(sj) * isn;
      rac3j[1] = cnr;
      rac3j[0] *= cnr;
      // returns
    } else { // label 20
      int nm = nmmo + 1;
      int nmat = (nm + 1) / 2;
      rac3j[nm - 1] = 1.0;
      rac3j[nmmo - 1] = dm / cm;
      double sjt = 1.0;
      double sjr = 1.0;
      if (nmat != nmmo) {
	int nbr = nmmo - nmat;
	for (int ibr45 = 1; ibr45 <= nbr; ibr45++) {
	  int irr = nm - ibr45;
	  m2--;
	  m3 = m1 + m2;
	  double cmp = cm;
	  cm = sqrt(1.0 * (j2po - m2) * (j2 + m2) * (j3po - m3) * (j3 + m3));
	  sjt = rac3j[irr - 1] * rac3j[irr - 1];
	  dm = 1.0 * (id + m2 * m3 * 2);
	  rac3j[irr - 1] *= ((dm - rac3j[irr] * cmp) / cm);
	  sjr += sjt;
	} // ibr45 loop
      }
      // label 50
      double osjt = sjt;
      sjt = rac3j[nmat - 1] * rac3j[nmat - 1];
      if (sjt >= osjt) {
	sjr += sjt;
      } else { // label 55
	nmat++;
      }
      // label 60
      double racmat = rac3j[nmat - 1];
      rac3j[0] = 1.0;
      m2 = mmn;
      m3 = m1 + m2;
      double cmp = sqrt(1.0 * (j2 - m2) * (j2po + m2) * (j3 - m3) * (j3po + m3));
      dm = 1.0 * (id + m2 * m3 * 2);
      rac3j[1] = dm / cmp;
      double sjl = 1.0;
      int nmatmo = nmat - 1;
      if (nmatmo > 1) {
	for (int irl70 = 2; irl70 <= nmatmo; irl70++) {
	  m2++;
	  m3 = m1 + m2;
	  cm = cmp;
	  cmp = sqrt(1.0 * (j2 - m2) * (j2po + m2) * (j3 - m3) * (j3po + m3));
	  sjt = rac3j[irl70 - 1] * rac3j[irl70 - 1];
	  dm = 1.0 * (id + m2 * m3 * 2);
	  rac3j[irl70] = (rac3j[irl70 - 1] * dm - rac3j[irl70 - 2] * cm) / cmp;
	  sjl += sjt;
	}
      } // label 75
      double ratrac = racmat / rac3j[nmat - 1];
      double rats = ratrac * ratrac;
      double sj = (sjr + sjl * rats) * j1tpo;
      rac3j[nmat - 1] = racmat;
      double cnr = 1.0 / sqrt(sj) * isn;
      for (int irr80 = nmat; irr80 <= nm; irr80++) rac3j[irr80 - 1] *= cnr;
      double cnl = cnr * ratrac;
      for (int irl85 = 1; irl85 <= nmatmo; irl85++) rac3j[irl85 - 1] *= cnl;
      // returns
    }
  }
}

void raba(
	  int le, dcomplex **am0m, dcomplex **w, double **tqce,
	  dcomplex **tqcpe, double **tqcs, dcomplex **tqcps
) {
  dcomplex **a, **ctqce, **ctqcs;
  // dcomplex c3;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  const double sq2i = 1.0 / sqrt(2.0);
  int nlem = le * (le + 2);
  const int nlemt = nlem + nlem;
  a = new dcomplex*[nlemt];
  ctqce = new dcomplex*[2];
  ctqcs = new dcomplex*[2];
  a[0] = new dcomplex[nlemt*2]();
  dcomplex *vec_a = a[0];
  ctqce[0] = new dcomplex[6]();
  ctqcs[0] = new dcomplex[6]();
  for (int ai = 1; ai < nlemt; ai++) a[ai] = a[0]+ai*2;
  ctqce[1] = ctqce[0]+3;
  ctqcs[1] = ctqcs[0]+3;
  dcomplex *vec_am0m = am0m[0];
  dcomplex *vec_w = w[0];
#ifdef USE_NVTX
  nvtxRangePush("raba outer loop 1");
#endif
#pragma omp parallel for
  for (int i20 = 1; i20 <= nlemt; i20++) {
    int i = i20 - 1;
    dcomplex c1 = cc0;
    dcomplex c2 = cc0;
#ifdef USE_NVTX
  nvtxRangePush("raba inner loop 1");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:c1, c2)
// #else
// #pragma omp parallel for simd reduction(+:c1, c2)
// #endif
#pragma omp parallel for simd reduction(+:c1, c2)
  for (int j10 = 1; j10 <= nlemt; j10++) {
      int j = j10 - 1;
      c1 += (vec_am0m[i*nlemt+j] * vec_w[4*j]);
      c2 += (vec_am0m[i*nlemt+j] * vec_w[4*j+1]);
    } // j10 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
    vec_a[2*i] = c1;
    vec_a[2*i+1] = c2;
  } //i20 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_NVTX
  nvtxRangePush("raba outer loop 2");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp teams distribute parallel for
// #else
// #pragma omp parallel for
// #endif
#pragma omp parallel for
  for (int ipo = 0; ipo < 2; ipo++) {
    int jpo = 1 - ipo;
    ctqce[ipo][0] = cc0;
    ctqce[ipo][1] = cc0;
    ctqce[ipo][2] = cc0;
    tqcpe[ipo][0] = cc0;
    tqcpe[ipo][1] = cc0;
    tqcpe[ipo][2] = cc0;
    ctqcs[ipo][0] = cc0;
    ctqcs[ipo][1] = cc0;
    ctqcs[ipo][2] = cc0;
    tqcps[ipo][0] = cc0;
    tqcps[ipo][1] = cc0;
    tqcps[ipo][2] = cc0;
    dcomplex &ctqce0 = ctqce[ipo][0];
    dcomplex &ctqce1 = ctqce[ipo][1];
    dcomplex &ctqce2 = ctqce[ipo][2];
    dcomplex &tqcpe0 = tqcpe[ipo][0];
    dcomplex &tqcpe1 = tqcpe[ipo][1];
    dcomplex &tqcpe2 = tqcpe[ipo][2];
    dcomplex &ctqcs0 = ctqcs[ipo][0];
    dcomplex &ctqcs1 = ctqcs[ipo][1];
    dcomplex &ctqcs2 = ctqcs[ipo][2];
    dcomplex &tqcps0 = tqcps[ipo][0];
    dcomplex &tqcps1 = tqcps[ipo][1];
    dcomplex &tqcps2 = tqcps[ipo][2];
    int kmax = le*(le+2);
    // for efficiency I should also linearise array w, but I cannot easily since I do not know for sure its major dimension (changes to containing class needed)
#ifdef USE_NVTX
    nvtxRangePush("raba inner loop 2");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:ctqce0, ctqce1, ctqce2, ctqcs0, ctqcs1, ctqcs2, tqcpe0, tqcpe1, tqcpe2, tqcps0, tqcps1, tqcps2)
// #else
// #pragma omp parallel for simd reduction(+:ctqce0, ctqce1, ctqce2, ctqcs0, ctqcs1, ctqcs2, tqcpe0, tqcpe1, tqcpe2, tqcps0, tqcps1, tqcps2)
// #endif
#pragma omp parallel for simd reduction(+:ctqce0, ctqce1, ctqce2, ctqcs0, ctqcs1, ctqcs2, tqcpe0, tqcpe1, tqcpe2, tqcps0, tqcps1, tqcps2)
    for (int k = 1; k<=kmax; k++) {
      int l60 = (int) sqrt(k+1);
      int im60 = k - (l60*l60) + 1;
      if (im60 == 0) {
	l60--;
	im60 = 2*l60+1;
      }
      else if (im60 > 2*l60 + 1) {
	im60 -= 2*l60 + 1;
	l60++;
      }
      int lpo = l60 + 1;
      int il = l60 * lpo;
      int m = im60 - lpo;
      int i = m + il - 1;
      int ie = i + nlem;
      int mmmu = m + 1;
      int mmmmu = (mmmu > 0) ? mmmu : -mmmu;
      double  rmu = 0.0;
      dcomplex acw;
      dcomplex acwp;
      dcomplex aca;
      dcomplex acap;
      if (mmmmu <= l60) {
	int immu = mmmu + il - 1;
	int immue = immu + nlem;
	rmu = -sqrt(1.0 * (l60 + mmmu) * (l60 - m)) * sq2i;
	acw = dconjg(vec_a[2*i+ipo]) * vec_w[4*immu+ipo] + dconjg(vec_a[2*ie+ipo]) * vec_w[4*immue+ipo];
	acwp = dconjg(vec_a[2*i+ipo]) * vec_w[4*immu+jpo] + dconjg(vec_a[2*ie+ipo]) * vec_w[4*immue+jpo];
	aca = dconjg(vec_a[2*i+ipo]) * vec_a[2*immu+ipo] + dconjg(vec_a[2*ie+ipo]) * vec_a[2*immue+ipo];
	acap = dconjg(vec_a[2*i+ipo]) * vec_a[2*immu+jpo] + dconjg(vec_a[2*ie+ipo]) * vec_a[2*immue+jpo];
	ctqce0 += (acw * rmu);
	tqcpe0 += (acwp * rmu);
	ctqcs0 += (aca * rmu);
	tqcps0 += (acap * rmu);
      }
      // label 30
      rmu = -1.0 * m;
      acw = dconjg(vec_a[2*i+ipo]) * vec_w[4*i+ipo] + dconjg(vec_a[2*ie+ipo]) * vec_w[4*ie+ipo];
      acwp = dconjg(vec_a[2*i+ipo]) * vec_w[4*i+jpo] + dconjg(vec_a[2*ie+ipo]) * vec_w[4*ie+jpo];
      aca = dconjg(vec_a[2*i+ipo]) * vec_a[2*i+ipo] + dconjg(vec_a[2*ie+ipo]) * vec_a[2*ie+ipo];
      acap = dconjg(vec_a[2*i+ipo]) * vec_a[2*i+jpo] + dconjg(vec_a[2*ie+ipo]) * vec_a[2*ie+jpo];
      ctqce1 += (acw * rmu);
      tqcpe1 += (acwp * rmu);
      ctqcs1 += (aca * rmu);
      tqcps1 += (acap * rmu);
      mmmu = m - 1;
      mmmmu = (mmmu > 0) ? mmmu : -mmmu;
      if (mmmmu <= l60) {
	int immu = mmmu + il - 1;
	int immue = immu + nlem;
	rmu = sqrt(1.0 * (l60 - mmmu) * (l60 + m)) * sq2i;
	acw = dconjg(vec_a[2*i+ipo]) * vec_w[4*immu+ipo] + dconjg(vec_a[2*ie+ipo]) * vec_w[4*immue+ipo];
	acwp = dconjg(vec_a[2*i+ipo]) * vec_w[4*immu+jpo] + dconjg(vec_a[2*ie+ipo]) * vec_w[4*immue+jpo];
	aca = dconjg(vec_a[2*i+ipo]) * vec_a[2*immu+ipo] + dconjg(vec_a[2*ie+ipo]) * vec_a[2*immue+ipo];
	acap = dconjg(vec_a[2*i+ipo]) * vec_a[2*immu+jpo] + dconjg(vec_a[2*ie+ipo]) * vec_a[2*immue+jpo];
	ctqce2 += (acw * rmu);
	tqcpe2 += (acwp * rmu);
	ctqcs2 += (aca * rmu);
	tqcps2 += (acap * rmu);
      } // ends if clause
    } // k loop (previously the l60 and im60 loops
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  } // ipo70 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_NVTX
  nvtxRangePush("raba loop 3");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd
// #else
// #pragma omp parallel for simd
// #endif
#pragma omp parallel for simd
  for (int ipo78 = 1; ipo78 <= 2; ipo78++) {
    int ipo = ipo78 - 1;
    tqce[ipo][0] = real(ctqce[ipo][0] - ctqce[ipo][2]) * sq2i;
    tqce[ipo][1] = real((ctqce[ipo][0] + ctqce[ipo][2]) * uim) * sq2i;
    tqce[ipo][2] = real(ctqce[ipo][1]);
    dcomplex c1 = tqcpe[ipo][0];
    dcomplex c2 = tqcpe[ipo][1];
    dcomplex c3 = tqcpe[ipo][2];
    tqcpe[ipo][0] = (c1 - c3) * sq2i;
    tqcpe[ipo][1] = (c1 + c3) * (uim * sq2i);
    tqcpe[ipo][2] = c2;
    tqcs[ipo][0] = -sq2i * real(ctqcs[ipo][0] - ctqcs[ipo][2]);
    tqcs[ipo][1] = -sq2i * real((ctqcs[ipo][0] + ctqcs[ipo][2]) * uim);
    tqcs[ipo][2] = -1.0 * real(ctqcs[ipo][1]);
    c1 = tqcps[ipo][0];
    c2 = tqcps[ipo][1];
    c3 = tqcps[ipo][2];
    tqcps[ipo][0] = -(c1 - c3) * sq2i;
    tqcps[ipo][1] = -(c1 + c3) * (uim * sq2i);
    tqcps[ipo][2] = -c2;
  } // ipo78 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  delete[] a[0];
  delete[] ctqce[0];
  delete[] ctqcs[0];
  delete[] a;
  delete[] ctqce;
  delete[] ctqcs;
}

void rftr(
	  double *u, double *up, double *un, double *gapv, double extins, double scatts,
	  double &rapr, double &cosav, double &fp, double &fn, double &fk, double &fx,
	  double &fy, double &fz
) {
  fk = u[0] * gapv[0] + u[1] * gapv[1] + u[2] * gapv[2];
  rapr = extins - fk;
  cosav = fk / scatts;
  fp = -(up[0] * gapv[0] + up[1] * gapv[1] + up[2] * gapv[2]);
  fn = -(un[0] * gapv[0] + un[1] * gapv[1] + un[2] * gapv[2]);
  fk = rapr;
  fx = u[0] * extins - gapv[0];
  fy = u[1] * extins - gapv[1];
  fz = u[2] * extins - gapv[2];
}

void scr0(double vk, double exri, ParticleDescriptor *c1) {
  const dcomplex cc0 = 0.0 + 0.0 * I;
  double exdc = exri * exri;
  double ccs = 4.0 * acos(0.0) / (vk * vk);
  double cccs = ccs / exdc;
  dcomplex csam = -(ccs / (exri * vk)) * 0.5 * I;
  //double scs = 0.0, ecs = 0.0, acs = 0.0;
  dcomplex *vec_rmi = c1->rmi[0];
  dcomplex *vec_rei = c1->rei[0];
#ifdef USE_NVTX
  nvtxRangePush("scr0 outer loop 1");
#endif
#pragma omp parallel for
  for (int i14 = 1; i14 <= c1->nsph; i14++) {
    int iogi = c1->iog[i14 - 1];
    if (iogi >= i14) {
      double sums = 0.0;
      dcomplex sum21 = cc0;
#ifdef USE_NVTX
      nvtxRangePush("scr0 inner loop 1");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:sums, sum21)
// #else
// #pragma omp parallel for simd reduction(+:sums, sum21)
// #endif
#pragma omp parallel for simd reduction(+:sums, sum21)
      for (int l10 = 1; l10 <= c1->li; l10++) {
	double fl = 1.0 * (l10 + l10 + 1);
	// dcomplex rm = 1.0 / c1->rmi[l10 - 1][i14 - 1];
	// dcomplex re = 1.0 / c1->rei[l10 - 1][i14 - 1];
	int vecindex = (l10 - 1)*c1->nsph + i14 - 1;
	dcomplex rm = 1.0 / vec_rmi[vecindex];
	dcomplex re = 1.0 / vec_rei[vecindex];
	double rvalue = real(dconjg(rm) * rm + dconjg(re) * re) * fl;
	sums += rvalue;
	sum21 += ((rm + re) * fl);
      } // l10 loop
#ifdef USE_NVTX
      nvtxRangePop();
#endif
      sum21 *= -1.0;
      double scasec = cccs * sums;
      double extsec = -cccs * real(sum21);
      double abssec = extsec - scasec;
      c1->sscs[i14 - 1] = scasec;
      c1->sexs[i14 - 1] = extsec;
      c1->sabs[i14 - 1] = abssec;
      double gcss = c1->gcsv[i14 - 1];
      c1->sqscs[i14 - 1] = scasec / gcss;
      c1->sqexs[i14 - 1] = extsec / gcss;
      c1->sqabs[i14 - 1] = abssec / gcss;
      c1->fsas[i14 - 1] = sum21 * csam;
    }
    // label 12
    // scs += c1->sscs[iogi - 1];
    // ecs += c1->sexs[iogi - 1];
    // acs += c1->sabs[iogi - 1];
    // tfsas += c1->fsas[iogi - 1];
  } // i14 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  double scs = 0.0;
  double ecs = 0.0;
  double acs = 0.0;
  dcomplex tfsas = cc0;
#ifdef USE_NVTX
  nvtxRangePush("scr0 loop 2");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:scs, ecs, acs, tfsas)
// #else
// #pragma omp parallel for simd reduction(+:scs, ecs, acs, tfsas)
// #endif
#pragma omp parallel for simd reduction(+:scs, ecs, acs, tfsas)
  for (int i14 = 1; i14 <= c1->nsph; i14++) {
    int iogi = c1->iog[i14 - 1];
    scs += c1->sscs[iogi - 1];
    ecs += c1->sexs[iogi - 1];
    acs += c1->sabs[iogi - 1];
    tfsas += c1->fsas[iogi - 1];
  }
  c1->scs = scs;
  c1->ecs = ecs;
  c1->acs = acs;
  c1->tfsas = tfsas;
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}

void scr2(
	  double vk, double vkarg, double exri, double *duk,
	  ParticleDescriptor *c1
) {
#ifdef USE_NVTX
  nvtxRangePush("scr2 starts");
#endif
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const dcomplex uim = 0.0 + 1.0 * I;
  dcomplex s11, s21, s12, s22, rm, re, csam, cph, phas, cc;
  double ccs = 1.0 / (vk * vk);
  csam = -(ccs / (exri * vk)) * 0.5 * I;
  const double pi4sq = 64.0 * acos(0.0) * acos(0.0);
  double cfsq = 4.0 / (pi4sq * ccs * ccs);
  cph = uim * exri * vkarg;
  int ls = (c1->li < c1->le) ? c1->li : c1->le;
  int kmax = ls*(ls+2);
  dcomplex *vec_rmi = c1->rmi[0];
  dcomplex *vec_rei = c1->rei[0];
  dcomplex *vec_w = c1->w[0];
  dcomplex *vec_sas = c1->sas[0][0];
#ifdef USE_NVTX
  nvtxRangePush("scr2 outer loop 1");
#endif
#pragma omp parallel for
  for (int i14 = 1; i14 <= c1->nsph; i14++) {
    int i = i14 - 1;
    int iogi = c1->iog[i14 - 1];
    if (iogi >= i14) {
      // int k = 0;
      dcomplex s11 = cc0;
      dcomplex s21 = cc0;
      dcomplex s12 = cc0;
      dcomplex s22 = cc0;
      // To parallelise, I run a linearised loop directly over k
      // working out the algebra, it turns out that
      // k = l10*l10-1+im10
      // we invert this to find
      // l10 = (int) sqrt(k+1) and im10 = k - l10*10+1
      // but if it results im10 = 0, then we set l10 = l10-1 and im10 = 2*l10+1
      // furthermore if it results im10 > 2*l10+1, then we set
      // im10 = im10 -(2*l10+1) and l10 = l10+1 (there was a rounding error in a nearly exact root)
#ifdef USE_NVTX
      nvtxRangePush("scr2 inner loop 1");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(-:s11, s21, s12, s22)
// #else
// #pragma omp parallel for simd reduction(-:s11, s21, s12, s22)
// #endif
#pragma omp parallel for simd reduction(-:s11, s21, s12, s22)
      for (int k = 1; k<=kmax; k++) {
	int l10 = (int) sqrt(k+1);
	int im10 = k - (l10*l10) + 1;
	if (im10 == 0) {
	  l10--;
	  im10 = 2*l10+1;
	}
	else if (im10 > 2*l10 + 1) {
	  im10 -= 2*l10 + 1;
	  l10++;
	}
	// I have all the indices in my linearised loop
	int l =  l10 - 1;
	int ke = k + c1->nlem;
	int km1t4 = (k - 1)*4;
	int kem1t4 = (ke - 1)*4;
	dcomplex auxrm0 = vec_w[km1t4] / vec_rmi[l*c1->nsph+i];
	dcomplex auxrm1 = vec_w[km1t4+1] / vec_rmi[l*c1->nsph+i];
	dcomplex auxre0 = vec_w[kem1t4] / vec_rei[l*c1->nsph+i];
	dcomplex auxre1 = vec_w[kem1t4+1] / vec_rei[l*c1->nsph+i];
	s11 -= vec_w[km1t4+2] * auxrm0 + vec_w[kem1t4+2] * auxre0;
	s21 -= vec_w[km1t4+3] * auxrm0 + vec_w[kem1t4+3] * auxre0;
	s12 -= vec_w[km1t4+2] * auxrm1 + vec_w[kem1t4+2] * auxre1;
	s22 -= vec_w[km1t4+3] * auxrm1 + vec_w[kem1t4+3] * auxre1;
      }
#ifdef USE_NVTX
      nvtxRangePop();
#endif
      int vecindex = i*4;
      vec_sas[vecindex] = s11 * csam;
      vec_sas[vecindex+2] = s21 * csam;
      vec_sas[vecindex+1] = s12 * csam;
      vec_sas[vecindex+3] = s22 * csam;
    }
    // label 12
    // dcomplex phas = cexp(cph * (duk[0] * c1->rxx[i] + duk[1] * c1->ryy[i] + duk[2] * c1->rzz[i]));
    // tsas00 += (c1->sas[iogi - 1][0][0] * phas);
    // tsas10 += (c1->sas[iogi - 1][1][0] * phas);
    // tsas01 += (c1->sas[iogi - 1][0][1] * phas);
    // tsas11 += (c1->sas[iogi - 1][1][1] * phas);
  } // i14 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  dcomplex tsas00 = cc0;
  dcomplex tsas10 = cc0;
  dcomplex tsas01 = cc0;
  dcomplex tsas11 = cc0;
#ifdef USE_NVTX
  nvtxRangePush("scr2 loop 2");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd reduction(+:tsas00, tsas10, tsas01, tsas11)
// #else
// #pragma omp parallel for simd reduction(+:tsas00, tsas10, tsas01, tsas11)
// #endif
#pragma omp parallel for simd reduction(+:tsas00, tsas10, tsas01, tsas11)
  for (int i14 = 1; i14 <= c1->nsph; i14++) {
    int i = i14 - 1;
    int iogi = c1->iog[i14 - 1];
    // label 12
    dcomplex phas = cexp(cph * (duk[0] * c1->rxx[i] + duk[1] * c1->ryy[i] + duk[2] * c1->rzz[i]));
    tsas00 += (c1->sas[iogi - 1][0][0] * phas);
    tsas10 += (c1->sas[iogi - 1][1][0] * phas);
    tsas01 += (c1->sas[iogi - 1][0][1] * phas);
    tsas11 += (c1->sas[iogi - 1][1][1] * phas);
  } // i14 loop
  c1->tsas[0][0] = tsas00;
  c1->tsas[1][0] = tsas10;
  c1->tsas[0][1] = tsas01;
  c1->tsas[1][1] = tsas11;
#ifdef USE_NVTX
  nvtxRangePop();
  //#endif
  //dcomplex *vec_vints = c1->vints[0];
  //#ifdef USE_NVTX
  nvtxRangePush("scr2 outer loop 3");
#endif
#pragma omp parallel for
  for (int i24 = 1; i24 <= c1->nsph; i24++) {
    int iogi = c1->iog[i24 - 1];
    if (iogi >= i24) {
      // int j = 0;
#ifdef USE_NVTX
      nvtxRangePush("scr2 inner loop 3");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd collapse(4)
// #else
// #pragma omp parallel for simd collapse(4)
// #endif
#pragma omp parallel for simd collapse(4)
      for (int ipo1 = 1; ipo1 <=2; ipo1++) {
	for (int jpo1 = 1; jpo1 <= 2; jpo1++) {
	  for (int ipo2 = 1; ipo2 <= 2; ipo2++) {
	    for (int jpo2 = 1; jpo2 <= 2; jpo2++) {
	      int j = jpo2 - 1 + (ipo2 - 1) * 2 + (jpo1 - 1) * 4 + (ipo1 - 1) * 8;
	      c1->vints[i24 - 1][j] = c1->sas[i24 - 1][jpo2 - 1][ipo2 - 1] * dconjg(c1->sas[i24 - 1][jpo1 - 1][ipo1 - 1]) * cfsq;
	    } // jpo2 loop
	  } // ipo2 loop
	} // jpo1 loop
      } // ipo1 loop
#ifdef USE_NVTX
      nvtxRangePop();
#endif
    }
  } // i24 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
  // int j = 0;
#ifdef USE_NVTX
  nvtxRangePush("scr2 loop 4");
#endif
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for collapse(4)
// #else
// #pragma omp parallel for collapse(4)
// #endif
#pragma omp parallel for collapse(4)
  for (int ipo1 = 1; ipo1 <=2; ipo1++) {
    for (int jpo1 = 1; jpo1 <= 2; jpo1++) {
      for (int ipo2 = 1; ipo2 <= 2; ipo2++) {
	for (int jpo2 = 1; jpo2 <= 2; jpo2++) {
	  int j = jpo2-1 + (ipo2-1)*2 + (jpo1-1)*4 + (ipo1-1)*8;
	  c1->vintt[j] = c1->tsas[jpo2 - 1][ipo2 - 1] * dconjg(c1->tsas[jpo1 - 1][ipo1 - 1]) * cfsq;
	} // jpo2 loop
      } // ipo2 loop
    } // jpo1 loop
  } // ipo1 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}

void str(double **rcf, ParticleDescriptor *c1) {
  dcomplex *ylm;
  const double pi = acos(-1.0);
  int last_configuration;
  c1->gcs = 0.0;
  double gcss = 0.0;
  last_configuration = 0;
  for (int i18 = 1; i18 <= c1->nsph; i18++) {
    int iogi = c1->iog[i18 - 1];
    if (iogi >= i18) {
      last_configuration++;
      gcss = pi * c1->ros[last_configuration - 1] * c1->ros[last_configuration - 1];
      c1->gcsv[i18 - 1] = gcss;
      int nsh = c1->nshl[last_configuration - 1];
      for (int j16 = 1; j16 <= nsh; j16++) {
	c1->rc[last_configuration - 1][j16 - 1] = rcf[last_configuration - 1][j16 - 1] * c1->ros[last_configuration - 1];
      } // j16 loop
    }
    c1->gcs += gcss;
  } // i18 loop
  int ylm_size = (c1->litpos > c1->lmtpos) ? c1->litpos : c1->lmtpos;
  ylm = new dcomplex[ylm_size]();
  int i = 0;
  for (int l1po28 = 1; l1po28 <= c1->lmpo; l1po28++) {
    int l1 = l1po28 - 1;
    for (int l2 = 1; l2 <= c1->lm; l2++) {
      r3j000(l1, l2, c1->rac3j);
      c1->ind3j[l1po28 - 1][l2 - 1] = i;
      int lmnpo = (l2 > l1) ? l2 - l1 + 1 : l1 - l2 + 1;
      int lmxpo = l2 + l1 + 1;
      int lpo28 = lmnpo;
      int il = 0;
      while (lpo28 <= lmxpo) {
	i++;
	il++;
	c1->v3j0[i - 1] = c1->rac3j[il - 1];
	lpo28 += 2;
      }
    } // l2 loop
  } // l1po28 loop
  int nsphmo = c1->nsph - 1;
  int lit = c1->li + c1->li;
  int ivy = 0;
  for (int nf40 = 1; nf40 <= nsphmo; nf40++) { // GPU portable?
    int nfpo = nf40 + 1;
    for (int ns40 = nfpo; ns40 <= c1->nsph; ns40++) {
      double rx = c1->rxx[nf40 - 1] - c1->rxx[ns40 - 1];
      double ry = c1->ryy[nf40 - 1] - c1->ryy[ns40 - 1];
      double rz = c1->rzz[nf40 - 1] - c1->rzz[ns40 - 1];
      double rr = 0.0;
      double crth = 0.0, srth = 0.0, crph = 0.0, srph = 0.0;
      polar(rx, ry, rz, rr, crth, srth, crph, srph);
      sphar(crth, srth, crph, srph, lit, ylm);
      for (int iv38 = 1; iv38 <= c1->litpos; iv38++) {
	c1->vyhj[iv38 + ivy - 1] = dconjg(ylm[iv38 - 1]);
      } // iv38 loop
      ivy += c1->litpos;
    } // ns40 loop
  } // nf40 loop
  int lmt = c1->li + c1->le;
  ivy = 0;
  for (int nf50 = 1; nf50 <= c1->nsph; nf50++) {
    double rx = c1->rxx[nf50 - 1];
    double ry = c1->ryy[nf50 - 1];
    double rz = c1->rzz[nf50 - 1];
    if (rx != 0.0 || ry != 0.0 || rz != 0.0) {
      double rr = 0.0;
      double crth = 0.0, srth = 0.0, crph = 0.0, srph = 0.0;
      polar(rx, ry, rz, rr, crth, srth, crph, srph);
      sphar(crth, srth, crph, srph, lmt, ylm);
      for (int iv48 = 1; iv48 <= c1->lmtpos; iv48++) {
	c1->vyj0[iv48 + ivy - 1] = dconjg(ylm[iv48 - 1]);
      } // iv48 loop
    }
    ivy += c1->lmtpos;
  } // nf50 loop
  delete[] ylm;
}

void tqr(
	 double *u, double *up, double *un, double *tqev, double *tqsv, double &tep,
	 double &ten, double &tek, double &tsp, double &tsn, double &tsk
) {
  tep = up[0] * tqev[0] + up[1] * tqev[1] + up[2] * tqev[2];
  ten = un[0] * tqev[0] + un[1] * tqev[1] + un[2] * tqev[2];
  tek = u[0] * tqev[0] + u[1] * tqev[1] + u[2] * tqev[2];
  tsp = up[0] * tqsv[0] + up[1] * tqsv[1] + up[2] * tqsv[2];
  tsn = un[0] * tqsv[0] + un[1] * tqsv[1] + un[2] * tqsv[2];
  tsk = u[0] * tqsv[0] + u[1] * tqsv[1] + u[2] * tqsv[2];
}

void ztm(dcomplex **am, ParticleDescriptor *c1) {
  // dcomplex gie, gle, a1, a2, a3, a4, sum1, sum2, sum3, sum4;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  // int i2 = 0; // old implementation
#ifdef USE_NVTX
  nvtxRangePush("ZTM starts");
#endif
#ifdef USE_NVTX
  nvtxRangePush("ZTM parallel loop 1");
#endif
  // C9 *c9_para = new C9(*c9);
  dcomplex *gis_v = c1->gis[0];
  dcomplex *gls_v = c1->gls[0];
  int k2max = c1->li*(c1->li+2);
  int k3max = c1->le*(c1->le+2);
  // To parallelise, I run a linearised loop directly over k
  // working out the algebra, it turns out that
  // k = l*l-1+im
  // we invert this to find
  // l = (int) sqrt(k+1) and im = k - l*l+1
  // but if it results im = 0, then we set l = l-1 and im = 2*l+1
  // furthermore if it results im > 2*l+1, then we set
  // im = im -(2*l+1) and l = l+1 (there was a rounding error in a nearly exact root)
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd collapse(3)
// #else
// #pragma omp parallel for simd collapse(3)
// #endif
#pragma omp parallel for simd collapse(3)
  for (int n2 = 1; n2 <= c1->nsph; n2++) { // GPU portable?
    for (int k2 = 1; k2<=k2max; k2++) {
      for (int k3 = 1; k3<=k3max; k3++) {
	int l2 = (int) sqrt(k2+1);
	int im2 = k2 - (l2*l2) + 1;
	if (im2 == 0) {
	  l2--;
	  im2 = 2*l2+1;
	}
	else if (im2 > 2*l2 + 1) {
	  im2 -= 2*l2 + 1;
	  l2++;
	}
	int l3 = (int) sqrt(k3+1);
	int im3 = k3 - (l3*l3) + 1;
	if (im3 == 0) {
	  l3--;
	  im3 = 2*l3+1;
	}
	else if (im3 > 2*l3 + 1) {
	  im3 -= 2*l3 + 1;
	  l3++;
	}
	// int l2tpo = l2 + l2 + 1;
	// int l3tpo = l3 + l3 + 1;
	int i2 = (n2-1) * c1->li * (c1->li + 2) + l2 * l2 + im2 - 1;
	int m2 = -l2 - 1 + im2;
	int i3 = l3 * l3 + im3 - 1;
	int m3 = -l3 - 1 + im3;
	int vecindex = (i2 - 1) * c1->nlem + i3 - 1;
	double *rac3j_local = (double *) malloc(c1->lmtpo*sizeof(double));
	gis_v[vecindex] = ghit_d(2, 0, n2, l2, m2, l3, m3, c1, rac3j_local);
	gls_v[vecindex] = ghit_d(2, 1, n2, l2, m2, l3, m3, c1, rac3j_local);
	free(rac3j_local);
      } // close k3 loop, former l3 + im3 loops
    } // close k2 loop, former l2 + im2 loops
  } // close n2 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_NVTX
  nvtxRangePush("ZTM loop 2");
#endif
  dcomplex *am_v = am[0];
  dcomplex *sam_v = c1->sam[0];
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd collapse(2)
// #else
// #pragma omp parallel for simd collapse(2)
// #endif
#pragma omp parallel for simd collapse(2)
  for (int i1 = 1; i1 <= c1->ndi; i1++) { // GPU portable?
    for (int i3 = 1; i3 <= c1->nlem; i3++) {
      dcomplex sum1 = cc0;
      dcomplex sum2 = cc0;
      dcomplex sum3 = cc0;
      dcomplex sum4 = cc0;
      int i1e = i1 + c1->ndi;
      int i3e = i3 + c1->nlem;
      for (int i2 = 1; i2 <= c1->ndi; i2++) {
	int i2e = i2 + c1->ndi;
	int vecindg_23 = (i2 - 1) * c1->nlem + i3 - 1;
	dcomplex gie = gis_v[vecindg_23];
	dcomplex gle = gls_v[vecindg_23];
	np_int vecinda_1 = (i1 - 1) * c1->ndit;
	np_int vecinda_1e = (i1 - 1 + c1->ndi) * c1->ndit;
	dcomplex a1 = am_v[vecinda_1 + i2 - 1];
	dcomplex a2 = am_v[vecinda_1 + i2e - 1];
	dcomplex a3 = am_v[vecinda_1e + i2 - 1];
	dcomplex a4 = am_v[vecinda_1e + i2e - 1];
	sum1 += (a1 * gie + a2 * gle);
	sum2 += (a1 * gle + a2 * gie);
	sum3 += (a3 * gie + a4 * gle);
	sum4 += (a3 * gle + a4 * gie);
      } // i2 loop
      int vecind1 = (i1 - 1) * c1->nlemt;
      int vecind1e = (i1e - 1) * c1->nlemt;
      sam_v[vecind1 + i3 - 1] = sum1;
      sam_v[vecind1 + i3e - 1] = sum2;
      sam_v[vecind1e + i3 - 1] = sum3;
      sam_v[vecind1e + i3e - 1] = sum4;
    } // i3 loop
  } // i1 loop
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd collapse(2)
// #else
// #pragma omp parallel for simd collapse(2)
// #endif 
#pragma omp parallel for simd collapse(2)
  for (int i1 = 1; i1 <= c1->ndi; i1++) {
    for (int i0 = 1; i0 <= c1->nlem; i0++) {
      int vecindex = (i1 - 1) * c1->nlem + i0 - 1;
      gis_v[vecindex] = dconjg(gis_v[vecindex]);
      gls_v[vecindex] = dconjg(gls_v[vecindex]);
    } // i0 loop
  } // i1 loop
  dcomplex *vec_am0m = c1->am0m[0];
// #ifdef USE_TARGET_OFFLOAD
// #pragma omp target teams distribute parallel for simd collapse(2)
// #else
// #pragma omp parallel for simd collapse(2)
// #endif
#pragma omp parallel for simd collapse(2)
  for (int i0 = 1; i0 <= c1->nlem; i0++) {
    for (int i3 = 1; i3 <= c1->nlemt; i3++) {
      int i0e = i0 + c1->nlem;
      dcomplex sum1 = cc0;
      dcomplex sum2 = cc0;
      for (int i1 = 1; i1 <= c1->ndi; i1 ++) {
	int i1e = i1 + c1->ndi;
	int vecind1 = (i1 - 1) * c1->nlemt;
	int vecind1e = (i1e - 1) * c1->nlemt;
	dcomplex a1 = sam_v[vecind1 + i3 - 1];
	dcomplex a2 = sam_v[vecind1e + i3 - 1];
	int vecindex = (i1 - 1) * c1->nlem + i0 - 1;
	dcomplex gie = gis_v[vecindex];
	dcomplex gle = gls_v[vecindex];
	sum1 += (a1 * gie + a2 * gle);
	sum2 += (a1 * gle + a2 * gie);
      } // i1 loop
      int vecind0 = (i0 - 1) * c1->nlemt;
      int vecind0e = (i0e - 1) * c1->nlemt;
      vec_am0m[vecind0 + i3 - 1] = -sum1;
      vec_am0m[vecind0e + i3 - 1] = -sum2;
      // c1->am0m[i0 - 1][i3 - 1] = -sum1;
      // c1->am0m[i0e - 1][i3 - 1] = -sum2;
    } // i3 loop
  } // i0 loop
#ifdef USE_NVTX
  nvtxRangePop();
#endif
#ifdef USE_NVTX
  nvtxRangePop();
#endif
}
