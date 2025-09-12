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

/*! \file inclu_subs.cpp
 *
 * \brief C++ implementation of INCLUSION subroutines.
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

using namespace std;

void cnf(int n, dcomplex z, int nm, dcomplex *csj, dcomplex *csy) {
  /*
  FROM CSPHJY OF LIBRARY specfun 
    
  ==========================================================
    Purpose: Compute spherical Bessel functions y
    Input :  Z --- Complex argument of y
             N --- Order of y ( N = 0,1,2,... )
             CSJ(N+1) --- j
             NM --- Highest order computed
    Output:  CSY(N+1) --- y
  ==========================================================
  */
  double a0 = cabs(z);
  if (a0 < 1.0e-60) {
    for (int k = 0; k <= n; k++) csy[k] = -1.0e300;
  } else {
    csy[0] = -ccos(z) / z;
    if (n > 0) {
      csy[1] = (csy[0] - csin(z)) / z;
      if (n > 1) {
	for (int k = 2; k <= nm; k++) {
	  double absjk = cabs(csj[k - 1]);
	  double absjkmo = cabs(csj[k - 2]);
	  csy[k] = (absjk > absjkmo) ?
	    (csj[k] * csy[k - 1] - 1.0 / (z * z)) / csj[k - 1] :
	    (csj[k] * csy[k - 2] - (2.0 * k - 1.0) / (z * z * z)) / csj[k - 2];
	}
      }
    }
  }
}

void exma(dcomplex **am, ParticleDescriptor *c1) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  dcomplex **at = c1->at;
  const int ndit = c1->ndit;
  const int ndm = c1->ndm;
  for (int j20 = 1; j20 <= c1->nlemt; j20++) {
    int j0 = j20 + ndit;
    for (int i20 = 1; i20 <= c1->nlemt; i20++) {
      dcomplex sum = cc0;
      for (int k = 0; k < ndm; k++) sum += at[i20 - 1][k] * am[k][j0 - 1];
      c1->am0m[i20 - 1][j20 - 1] = sum;
    }
  } // j20 loop
}

void incms(dcomplex **am, double enti, ParticleDescriptor *c1) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  dcomplex **at = c1->at;
  int nbl, i1;
  const int ndi = c1->ndi;
  const int ndit = ndi + ndi;
  const int ndm = c1->ndm;
  nbl = 0;
  for (int n1 = 1; n1 < c1->nsph; n1++) {
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
	  i1 = in1 + ilm1;
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
	      dcomplex cgh = ghit(0, 0, nbl, l1, m1, l2, m2, c1);
	      dcomplex cgk = ghit(0, 1, nbl, l1, m1, l2, m2, c1);
	      am[i1 - 1][i2 - 1] = cgh;
	      am[i1 - 1][i2e - 1] = cgk;
	      am[i1e - 1][i2 - 1] = cgk;
	      am[i1e - 1][i2e - 1] = cgh;
	      am[j1 - 1][j2 - 1] = cgh * ish;
	      am[j1 - 1][j2e - 1] = cgk * isk;
	      am[j1e - 1][j2 - 1] = cgk * isk;
	      am[j1e - 1][j2e - 1] = cgh * ish;
	    } // im2 loop 24
	  } // l2 loop 24
	} // im1 loop 24
      } // l1 loop 24
    } // n2 loop 26
  } // n1 loop 26
  for (int n1 = 1; n1 <= c1->nsph; n1++) {
    int in1 = (n1 - 1) * c1->nlim;
    for (int l1 = 1; l1 <= c1->li; l1++) {
      dcomplex frm = c1->rmi[l1 - 1][n1 - 1];
      dcomplex fre = c1->rei[l1 - 1][n1 - 1];
      int l1po = l1 + 1;
      int il1 = l1po * l1;
      int l1tpo = l1po + l1;
      for (int im1 = 1; im1 <= l1tpo; im1++) {
	int m1 = im1 - l1po;
	int ilm1 = il1 + m1;
	i1 = in1 + ilm1;
	int i1e = i1 + ndi;
	for (int ilm2 = 1; ilm2 <= c1->nlim; ilm2++) {
	  int i2 = in1 + ilm2;
	  int i2e = i2 + ndi;
	  am[i1 - 1][i2 - 1] = cc0;
	  am[i1 - 1][i2e - 1] = cc0;
	  am[i1e - 1][i2 - 1] = cc0;
	  am[i1e - 1][i2e - 1] = cc0;
	} // ilm2 loop 28
	am[i1 - 1][i1 - 1] = frm;
	am[i1e - 1][i1e - 1]= fre;
      } // im1 loop 30
    } // l1 loop 30
  } // n1 loop 30
  int nditpo = ndit + 1;
  for (i1 = 1; i1 <= c1->nlemt; i1++) {
    int i3 = i1 + ndit;
    for (int i2 = nditpo; i2 <= ndm; i2++) {
      am[i3 - 1][i2 - 1] = cc0;
      at[i1 - 1][i2 - 1] = cc0;
    } // i2 loop 40
  } // i1 loop 40
  i1 = 0;
  for (int l1 = 1; l1 <= c1->le; l1++) {
    dcomplex frm = c1->rm0[l1 - 1];
    dcomplex fre = c1->re0[l1 - 1];
    dcomplex ftm = c1->tm0[l1 - 1];
    dcomplex fte = c1->te0[l1 - 1];
    int l1tpo = l1 + l1 + 1;
    for (int im1  = 1; im1 <= l1tpo; im1 ++) {
      i1++;
      int i1e = i1 + c1->nlem;
      int i3 = i1 + ndit;
      int i3e = i3 + c1->nlem;
      am[i3 - 1][i3 - 1] = frm;
      am[i3e - 1][i3e - 1] = fre;
      at[i1 - 1][i3 - 1] = ftm;
      at[i1e - 1][i3e - 1] = fte;
    } // im1 loop 45
  } // l1 loop 45
  if (enti != 0.0) {
    for (int l2 = 1; l2 <= c1->le; l2++) {
      dcomplex frm = c1->rmw[l2 - 1];
      dcomplex fre = c1->rew[l2 - 1];
      dcomplex ftm = c1->tm[l2 - 1];
      dcomplex fte = c1->te[l2 - 1];
      int l2po = l2 + 1;
      int il2 = l2po * l2;
      int l2tpo = l2po + l2;
      for (int im2 = 1; im2 <= l2tpo; im2++) {
	int m2 = im2 - l2po;
	int i2 = il2 + m2;
	int j2 = il2 - m2;
	int i2e = i2 + c1->nlem;
	int j2e = j2 + c1->nlem;
	int i3 = i2 + ndit;
	int j3 = j2 + ndit;
	int i3e = i3 + c1->nlem;
	int j3e = j3 + c1->nlem;
	for (int n1 = 1; n1 <= c1->nsph; n1++) {
	  int in1 = (n1 - 1) * c1->nlim;
	  for (int l1 = 1; l1 <= c1->li; l1 ++) {
	    int l1po = l1 + 1;
	    int il1 = l1po * l1;
	    int l1tpo = l1po + l1;
	    for (int im1 = 1; im1 <= l1tpo; im1++) {
	      int m1 = im1 - l1po;
	      int ilm1 = il1 + m1;
	      int jlm1 = il1 - m1;
	      i1 = in1 + ilm1;
	      int i1e = i1 + ndi;
	      int j1 = in1 + jlm1;
	      int j1e = j1 + ndi;
	      int isil = ((m2 + m1) % 2 == 0) ? 1 : -1;
	      dcomplex cgi = ghit(2, 0, n1, l1, m1, l2, m2, c1);
	      dcomplex cgl = ghit(2, 1, n1, l1, m1, l2, m2, c1);
	      am[i1 - 1][i3 - 1] = cgi;
	      am[i1 - 1][i3e - 1] = cgl;
	      am[i1e - 1][i3 - 1] = cgl;
	      am[i1e - 1][i3e - 1] = cgi;
	      am[j3 - 1][j1 - 1] = cgi * frm * isil;
	      am[j3 - 1][j1e - 1] = cgl * frm * isil;
	      am[j3e - 1][j1 - 1] = cgl * fre * isil;
	      am[j3e - 1][j1e - 1] = cgi * fre * isil;
	      at[j2 - 1][j1 - 1] = cgi * ftm * isil;
	      at[j2 - 1][j1e - 1] = cgl * ftm * isil;
	      at[j2e - 1][j1 - 1] = cgl * fte * isil;
	      at[j2e - 1][j1e - 1] = cgi * fte * isil;
	      // returns
	    } // im1 loop 50
	  } // l1 loop 50
	} // n1 loop 50
      } // im2 loop 50
    } // l2 loop 50
  } else {
    // label 55
    int i2 = 0;
    for (int l2 = 1; l2 <= c1->le; l2++) {
      dcomplex frm = c1->rmw[l2 - 1];
      dcomplex fre = c1->rew[l2 - 1];
      dcomplex ftm = c1->tm[l2 - 1];
      dcomplex fte = c1->te[l2 - 1];
      int l2tpo = l2 + l2 + 1;
      int m2 = -l2 - 1;
      for (int im2 = 1; im2 <= l2tpo; im2++) {
	m2++;
	i2++;
	int i2e = i2 + c1->nlem;
	int i3 = i2 + ndit;
	int i3e = i3 + c1->nlem;
	i1 = 0;
	for (int n1 = 1; n1 <= c1->nsph; n1++) {
	  for (int l1 = 1; l1 <= c1->li; l1++) {
	    int l1tpo = l1 + l1 + 1;
	    int m1 = -l1 - 1;
	    for (int im1 = 1; im1 <= l1tpo; im1++) {
	      m1++;
	      i1++;
	      int i1e = i1 + ndi;
	      dcomplex cgi = ghit(2, 0, n1, l1, m1, l2, m2, c1);
	      dcomplex cgl = ghit(2, 1, n1, l1, m1, l2, m2, c1);
	      am[i1 - 1][i3 - 1] = cgi;
	      am[i1 - 1][i3e - 1] = cgl;
	      am[i1e - 1][i3 - 1] = cgl;
	      am[i1e - 1][i3e - 1] = cgi;
	      cgi = dconjg(cgi);
	      cgl = dconjg(cgl);
	      am[i3 - 1][i1 - 1] = cgi * frm;
	      am[i3 - 1][i1e - 1] = cgl * frm;
	      am[i3e - 1][i1 - 1] = cgl * fre;
	      am[i3e - 1][i1e - 1] = cgi * fre;
	      at[i2 - 1][i1 - 1] = cgi * ftm;
	      at[i2 - 1][i1e - 1] = cgl * ftm;
	      at[i2e - 1][i1 - 1] = cgl * fte;
	      at[i2e - 1][i1e - 1] = cgi * fte;
	    } // im1 loop 60
	  } // l1 loop 60
	} // n1 loop 60
      } // im2 loop 60
    } // l2 loop 60
  } // END OF enti = 0.0 CASE
}

void indme(
	   int i, int npnt, int npntts, double vk, dcomplex ent, double enti,
	   dcomplex entn, int &jer, int &lcalc, dcomplex &arg, ParticleDescriptor *c1) {
  const dcomplex uim = 0.0 + I * 1.0;
  const int nstp = npnt - 1;
  const int nstpts = npntts - 1;
  const int lipo = c1->li + 1;
  const int lipt = c1->li + 2;
  dcomplex *cfj = new dcomplex[lipt]();
  dcomplex *cfn = new dcomplex[lipt]();
  dcomplex *fb = new dcomplex[lipt]();
  dcomplex *fbi = new dcomplex[lipt]();
  dcomplex *fn = new dcomplex[lipt]();
  double *rfj = new double[lipt]();
  double *rfn = new double[lipt]();
  jer = 0;
  double sz = vk * c1->ros[i - 1];
  c1->vsz[i - 1] = sz;
  double vkr1 = vk * c1->rc[i - 1][0];
  int nsh = c1->nshl[i - 1];
  c1->vkt[i - 1] = csqrt(c1->dc0[0]);
  arg = vkr1 * c1->vkt[i - 1];
  dcomplex arin = arg;
  if (imag(arg) != 0.0) {
    cbf(lipo, arg, lcalc, cfj);
    if (lcalc < lipo) {
      jer = 5;
      delete[] cfj;
      delete[] cfn;
      delete[] fb;
      delete[] fbi;
      delete[] fn;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    // label 122
    for (int j124 = 0; j124 < lipt; j124++) fbi[j124] = cfj[j124];
  } else { // label 126
    double rarg = real(arg);
    rbf(lipo, rarg, lcalc, rfj);
    if (lcalc < lipo) {
      jer = 5;
      delete[] cfj;
      delete[] cfn;
      delete[] fb;
      delete[] fbi;
      delete[] fn;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    // label 128
    for (int j130 = 0; j130 < lipt; j130++) fbi[j130] = rfj[j130];
  }
  // label 132
  dcomplex aris = sz * entn;
  arg = aris;
  if (enti != 0.0) {
    cbf(lipo, arg, lcalc, cfj);
    if (lcalc < lipo) {
      jer = 11;
      delete[] cfj;
      delete[] cfn;
      delete[] fb;
      delete[] fbi;
      delete[] fn;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    cnf(lipo, arg, lcalc, cfj, cfn);
    // QUESTION: should we check for lcalc and throw JER=12 if failing?
    // see lines 2492 -2505 in INCLU.F (test done in REAL case but not in COMPLEX case).
    for (int j143 = 0; j143 < lipt; j143++) {
      fb[j143] = cfj[j143];
      fn[j143] = cfn[j143];
    }
  } else { // label 145
    double rarg = real(aris);
    rbf(lipo, rarg, lcalc, rfj);
    if (lcalc < lipo) {
      jer = 11;
      delete[] cfj;
      delete[] cfn;
      delete[] fb;
      delete[] fbi;
      delete[] fn;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    rnf(lipo, rarg, lcalc, rfn);
    if (lcalc < lipo) {
      jer = 12;
      delete[] cfj;
      delete[] cfn;
      delete[] fb;
      delete[] fbi;
      delete[] fn;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    for (int j150 = 0; j150 < lipt; j150++) {
      fb[j150] = rfj[j150];
      fn[j150] = rfn[j150];
    }
  }
  // label 152
  dcomplex *rmf = new dcomplex[c1->li]();
  dcomplex *drmf = new dcomplex[c1->li]();
  dcomplex *ref = new dcomplex[c1->li]();
  dcomplex *dref = new dcomplex[c1->li]();
  int ic = 0;
  if (nsh < 2) { // nsh == 1
    dcomplex cri = c1->dc0[0] / ent;
    for (int l160 = 1; l160 <= c1->li; l160++) {
      int lpo = l160 + 1;
      int ltpo = lpo + l160;
      int lpt = lpo + 1;
      dcomplex dfbi = (l160 * fbi[l160 - 1] - lpo * fbi[lpt - 1]) * arin + fbi[lpo - 1] * ltpo;
      dcomplex dfb = (l160 * fb[l160 - 1] - lpo * fb[lpt - 1]) * aris + fb[lpo - 1] * ltpo;
      dcomplex dfn = (l160 * fn[l160 - 1] - lpo * fn[lpt - 1]) * aris + fn[lpo - 1] * ltpo;
      dcomplex ccna = fbi[lpo - 1] * dfn;
      dcomplex ccnb = fn[lpo - 1] * dfbi;
      dcomplex ccnc = fbi[lpo - 1] * dfb;
      dcomplex ccnd = fb[lpo - 1] * dfbi;
      c1->rmi[l160 - 1][i - 1] = 1.0 + uim * (ccna - ccnb) / (ccnc - ccnd);
      c1->rei[l160 - 1][i - 1] = 1.0 + uim * (cri * ccna - ccnb) / (cri * ccnc - ccnd);
    } // l160 loop
  } else { // label 165: nsh > 1
    for (int l180 = 1; l180 <= c1->li; l180++) {
      int lpo = l180 + 1;
      int ltpo = lpo + l180;
      int lpt = lpo + 1;
      double dltpo = 1.0 * ltpo;
      dcomplex y1 = fbi[lpo - 1];
      dcomplex dy1 = (l180 * fbi[l180 - 1] - lpo * fbi[lpt - 1]) * c1->vkt[i - 1] / dltpo;
      dcomplex y2 = y1;
      dcomplex dy2 = dy1;
      ic = 0;
      for (int ns = 2; ns <= nsh; ns++) {
	int nsmo = ns - 1;
	double vkr = vk * c1->rc[i - 1][nsmo - 1];
	if (ns % 2 != 0) {
	  // ic is incremented before being read in this loop.
	  int step = vk * (c1->rc[i - 1][ns - 1] - c1->rc[i - 1][nsmo - 1]) / nstp;
	  arg = c1->dc0[++ic];
	  rkc(nstp, step, arg, vkr, lpo, y1, y2, dy1, dy2);
	} else { // label 170
	  diel(nstpts, nsmo, i, ic, vk, c1);
	  double stepts = vk * (c1->rc[i - 1][ns - 1] - c1->rc[i - 1][nsmo - 1]) / nstpts;
	  rkt(nstpts, stepts, vkr, lpo, y1, y2, dy1, dy2, c1);
	}
      } // ns loop 176
      rmf[l180 - 1] = y1 * sz;
      drmf[l180 - 1] = dy1 * sz + y1;
      ref[l180 - 1] = y2 * sz;
      dref[l180 - 1] = dy2 * sz + y2;
    } // l180 loop
    dcomplex cri = (nsh % 2 == 0) ? 1.0 + I * 0.0 : c1->dc0[ic - 1] / ent;
    for (int l190 = 1; l190 <= c1->li; l190++) {
      int lpo = l190 + 1;
      int ltpo = lpo + l190;
      int lpt = lpo + 1;
      dcomplex dfb = (l190 * fb[l190 - 1] - lpo * fb[lpt - 1]) * aris + fb[lpo - 1] * ltpo;
      dcomplex dfn = (l190 * fn[l190 - 1] - lpo * fn[lpt - 1]) * aris + fn[lpo - 1] * ltpo;
      dcomplex ccna = rmf[l190 - 1] * dfn;
      dcomplex ccnb = drmf[l190 - 1] * fn[lpo - 1] * sz * ltpo;
      dcomplex ccnc = rmf[l190 - 1] * dfb;
      dcomplex ccnd = drmf[l190 - 1] * fb[lpo - 1]* sz * ltpo;
      c1->rmi[l190 - 1][i -1] = 1.0 + uim * (ccna - ccnb) / (ccnc - ccnd);
      ccna = ref[l190 - 1] * dfn;
      ccnb = dref[l190 - 1] * fn[lpo] * sz * ltpo;
      ccnc = ref[l190 - 1] * dfb;
      ccnd = dref[l190 - 1] * fb[lpo - 1] * sz * ltpo;
      c1->rei[l190 - 1][i - 1] =1.0 + uim * (cri * ccna - ccnb) / (cri * ccnc - ccnd);
    } // l190 loop
  } // nsh if
  delete[] cfj;
  delete[] cfn;
  delete[] fb;
  delete[] fbi;
  delete[] fn;
  delete[] rfj;
  delete[] rfn;
  delete[] rmf;
  delete[] drmf;
  delete[] ref;
  delete[] dref;
}

void instr(double **rcf, ParticleDescriptor *c1) {
  const int ylm_size = (c1->litpos > c1->lmtpos) ? c1->litpos : c1->lmtpos;
  dcomplex *ylm = new dcomplex[ylm_size]();
  double rx, ry, rz, rr, crth, srth, crph, srph;
  int ivy;
  for (int i18 = 0; i18 < c1->nsph; i18++) {
    int i = i18 + 1;
    if (c1->iog[i18] >= i) {
      int nsh = c1->nshl[i18];
      for (int j = 0; j < nsh; j++)
	c1->rc[i18][j] = rcf[i18][j] * c1->ros[i18];
    }
  } // i18 loop
  int i = 0;
  for (int l1po = 1; l1po <= c1->lmpo; l1po++) {
    int l1 = l1po - 1;
    for (int l2 = 1; l2 <= c1->lm; l2++) {
      r3j000(l1, l2, c1->rac3j);
      c1->ind3j[l1po - 1][l2 - 1] = i;
      int lmnpo = 1 + ((l2 - l1 > 0) ? l2 - l1 : l1 - l2);
      int lmxpo = l2 + l1 + 1;
      int il = 0;
      int lpo = lmnpo;
      while (lpo <= lmxpo) {
	c1->v3j0[i++] = c1->rac3j[il++];
	lpo += 2;
      }
    } // l2 loop
  } // l1po loop 28
  int nsphmo = c1->nsph - 1;
  int lit = c1->li + c1->li;
  ivy = 0;
  for (int nf40 = 0; nf40 < nsphmo; nf40++) {
    int nf = nf40 + 1;
    for (int ns = nf; ns < c1->nsph; ns++) {
      rx = c1->rxx[nf40] - c1->rxx[ns];
      ry = c1->ryy[nf40] - c1->ryy[ns];
      rz = c1->rzz[nf40] - c1->rzz[ns];
      polar(rx, ry, rz, rr, crth, srth, crph, srph);
      sphar(crth, srth, crph, srph, lit, ylm);
      for (int iv38 = 0; iv38 < c1->litpos; iv38++)
	c1->vyhj[iv38 + ivy] = dconjg(ylm[iv38]);
      ivy += c1->litpos;
    } // ns loop
  } // nf40 loop
  int lmt = c1->li + c1->le;
  ivy = 0;
  for (int nf50 = 0; nf50 < c1->nsph; nf50++) {
    rx = c1->rxx[nf50];
    ry = c1->ryy[nf50];
    rz = c1->rzz[nf50];
    if (rx != 0.0 || ry != 0.0 || rz != 0.0) {
      polar(rx, ry, rz, rr, crth, srth, crph, srph);
      sphar(crth, srth, crph, srph, lmt, ylm);
      for (int iv48 = 0; iv48 < c1->lmtpos; iv48++)
	c1->vyj0[iv48 + ivy] = dconjg(ylm[iv48]);
      ivy += c1->lmtpos;
    }
  } // nf50 loop
  delete[] ylm;
}

void ospv(ParticleDescriptor *c1, double vk, double sze, double exri, dcomplex entn, double enti, int &jer, int &lcalc, dcomplex &arg) {
  const dcomplex uim = 0.0 + I * 1.0;
  const int nsph = c1->nsph;
  const int nsphmo = c1->nsph - 1;
  const int lit = c1->li + c1->li;
  const int litpo = lit + 1;
  const int array_size = (c1->litpo > c1->lmtpo) ? c1->litpo : c1->lmtpo;
  dcomplex *cfj = new dcomplex[array_size]();
  dcomplex *cfn = new dcomplex[array_size]();
  double *rfj = new double[array_size]();
  double *rfn = new double[array_size]();
  
  int ivhb = 0;
  for (int i130 = 1; i130 <= nsphmo; i130++) {
    int ipo = i130 + 1;
    for (int j130 = ipo; j130 <= nsph; j130++) {
      double rx = c1->rxx[j130 - 1] - c1->rxx[i130 - 1];
      double ry = c1->ryy[j130 - 1] - c1->ryy[i130 - 1];
      double rz = c1->rzz[j130 - 1] - c1->rzz[i130 - 1];
      double rr = sqrt(rx * rx + ry * ry + rz * rz);
      arg = rr * vk * entn;
      if (enti != 0.0) {
	cbf(lit, arg, lcalc, cfj);
	if (lcalc < lit) {
	  jer = 11;
	  delete[] cfj;
	  delete[] cfn;
	  delete[] rfj;
	  delete[] rfn;
	  return;
	}
	cnf(lit, arg, lcalc, cfj, cfn);
	for (int lpo = 0; lpo < litpo; lpo++) c1->vh[lpo + ivhb] = cfj[lpo] + uim * cfn[lpo];
	// goes to 130
      } else { // label 123
	double rarg = real(arg);
	rbf(lit, rarg, lcalc, rfj);
	if (lcalc < lit) {
	  jer = 11;
	  delete[] cfj;
	  delete[] cfn;
	  delete[] rfj;
	  delete[] rfn;
	  return;
	}
	rnf(lit, rarg, lcalc, rfn);
	if (lcalc < lit) {
	  jer = 12;
	  delete[] cfj;
	  delete[] cfn;
	  delete[] rfj;
	  delete[] rfn;
	  return;
	}
	for (int lpo = 0; lpo < litpo; lpo++) c1->vh[lpo + ivhb] = rfj[lpo] + uim * rfn[lpo];
      }
      // label 130
      ivhb += litpo;
    } // j130 loop
  } // i130 loop
  const int lmt = c1->li + c1->le;
  const int lmtpo = lmt + 1;
  ivhb = 0;
  for (int i155 = 1; i155 <= nsph; i155++) {
    double rx = c1->rxx[i155 - 1];
    double ry = c1->ryy[i155 - 1];
    double rz = c1->rzz[i155 - 1];
    if (rx != 0.0 || ry != 0.0 || rz != 0.0) {
      double rr = sqrt(rx * rx + ry * ry + rz * rz);
      arg = rr * vk * entn;
      if (enti != 0.0) {
	cbf(lmt, arg, lcalc, cfj);
	if (lcalc < lmt) {
	  jer = 11;
	  delete[] cfj;
	  delete[] cfn;
	  delete[] rfj;
	  delete[] rfn;
	  return;
	}
	for (int lpo = 0; lpo < lmtpo; lpo++) c1->vj0[lpo + ivhb] = cfj[lpo];
	// goes to 155
      } else { // label 150
	double rarg = real(arg);
	rbf(lmt, rarg, lcalc, rfj);
	if (lcalc < lmt) {
	  jer = 11;
	  delete[] cfj;
	  delete[] cfn;
	  delete[] rfj;
	  delete[] rfn;
	  return;
	}
	for (int lpo = 0; lpo < lmtpo; lpo++) c1->vj0[lpo + ivhb] = rfj[lpo];
      }
    }
    // label 155
    ivhb += lmtpo;
  } // i155 loop
  
  const int lepo = c1->le + 1;
  const int lept = c1->le + 2;
  dcomplex *fb0 = new dcomplex[lept]();
  dcomplex *fh0 = new dcomplex[lept]();
  dcomplex aris0 = sze * entn;
  arg = aris0;
  if (enti != 0.0) {
    cbf(lepo, arg, lcalc, cfj);
    if (lcalc < lepo) {
      jer = 11;
      delete[] cfj;
      delete[] cfn;
      delete[] rfj;
      delete[] rfn;
      delete[] fb0;
      delete[] fh0;
      return;
    }
    cnf(lepo, arg, lcalc, cfj, cfn);
    for (int j162 = 0; j162 < lept; j162++) {
      fb0[j162] = cfj[j162];
      fh0[j162] = cfj[j162] + uim * cfn[j162];
    } // j162 loop
    // goes to 170
  } else { // label 163
    double rarg = real(arg);
    rbf(lepo, rarg, lcalc, rfj);
    if (lcalc < lepo) {
      jer = 11;
      delete[] cfj;
      delete[] cfn;
      delete[] rfj;
      delete[] rfn;
      delete[] fb0;
      delete[] fh0;
      return;
    }
    rnf(lepo, rarg, lcalc, rfn);
    if (lcalc < lepo) {
      jer = 12;
      delete[] cfj;
      delete[] cfn;
      delete[] rfj;
      delete[] rfn;
      delete[] fb0;
      delete[] fh0;
      return;
    }
    for (int j168 = 0; j168 < lept; j168++) {
      fb0[j168] = rfj[j168];
      fh0[j168] = rfj[j168] + uim * rfn[j168];
    } // j168 loop
  }
  // label 170
  double arex = sze * exri;
  arg = arex;
  double rarg = arex;
  rbf(lepo, rarg, lcalc, rfj);
  if (lcalc < lepo) {
    jer = 1;
    delete[] cfj;
    delete[] cfn;
    delete[] rfj;
    delete[] rfn;
    delete[] fb0;
    delete[] fh0;
    return;
  }
  rnf(lepo, rarg, lcalc, rfn);
  if (lcalc < lepo) {
    jer = 2;
    delete[] cfj;
    delete[] cfn;
    delete[] rfj;
    delete[] rfn;
    delete[] fb0;
    delete[] fh0;
    return;
  }
  dcomplex *fbe = new dcomplex[lept]();
  dcomplex *fhe = new dcomplex[lept]();
  for (int j175 = 0; j175 < lept; j175++) {
    fbe[j175] = rfj[j175];
    fhe[j175] = rfj[j175] + uim * rfn[j175];
  } // j175 loop
  dcomplex cri = exri / entn;
  for (int l184 = 1; l184 <= c1->le; l184++) {
    int lpo = l184 + 1;
    int lpt = lpo + 1;
    double dltpo = 1.0 / (lpo + l184);
    dcomplex dfb0 = fb0[lpo - 1] + (l184 * fb0[l184 - 1] - lpo * fb0[lpt - 1]) * aris0 * dltpo;
    dcomplex dfh0 = fh0[lpo - 1] + (l184 * fh0[l184 - 1] - lpo * fh0[lpt - 1]) * aris0 * dltpo;
    dcomplex dfbe = fbe[lpo - 1] + (l184 * fbe[l184 - 1] - lpo * fbe[lpt - 1]) * arex * dltpo;
    dcomplex dfhe = fhe[lpo - 1] + (l184 * fhe[l184 - 1] - lpo * fhe[lpt - 1]) * arex * dltpo;
    dcomplex ccna = aris0 * fb0[lpo - 1] * dfhe;
    dcomplex ccnb = arex * fhe[lpo - 1] * dfb0;
    c1->rm0[l184 - 1] = -(ccna * cri - ccnb) * uim;
    c1->re0[l184 - 1] = -(ccna - ccnb * cri) * uim;
    ccna = aris0 * fh0[lpo - 1] * dfhe;
    ccnb = arex * fhe[lpo - 1] * dfh0;
    c1->rmw[l184 - 1] = -(ccna * cri - ccnb) * uim;
    c1->rew[l184 - 1] = -(ccna - ccnb * cri) * uim;
    ccna = aris0 * fh0[lpo - 1] * dfbe;
    ccnb = arex * fbe[lpo - 1] * dfh0;
    c1->tm[l184 - 1] = (ccna * cri - ccnb) * uim;
    c1->te[l184 - 1] = (ccna - ccnb * cri) * uim;
    ccna = aris0 * fb0[lpo - 1] * dfbe;
    ccnb = arex * fbe[lpo - 1] * dfb0;
    c1->tm0[l184 - 1] = (ccna * cri - ccnb) * uim;
    c1->te0[l184 - 1] = (ccna - ccnb * cri) * uim;
  } // l184 loop
  
  // Clean up memory.
  delete[] cfj;
  delete[] cfn;
  delete[] rfj;
  delete[] rfn;
  delete[] fb0;
  delete[] fh0;
  delete[] fbe;
  delete[] fhe;
}
