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

/*! \file sph_subs.cpp
 *
 * \brief C++ implementation of SPHERE subroutines.
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

void aps(double ****zpv, int li, int nsph, ParticleDescriptor *c1, double sqk, double *gaps) {
  dcomplex cc0 = 0.0 + 0.0 * I;
  dcomplex summ, sume, suem, suee, sum;
  double half_pi = acos(0.0);
  double cofs = half_pi * 2.0 / sqk;
  for (int i40 = 0; i40 < nsph; i40++) {
    int i = i40 + 1;
    int iogi = c1->iog[i40];
    if (iogi >= i) {
      sum = cc0;
      for (int l30 = 0; l30 < li; l30++) {
	int l = l30 + 1;
	int ltpo = l + l + 1;
	for (int ilmp = 1; ilmp < 4; ilmp++) {
	  int ilmp30 = ilmp - 1;
	  bool goto30 = (l == 1 && ilmp == 1) || (l == li && ilmp == 3);
	  if (!goto30) {
	    int lmpml = ilmp - 2;
	    int lmp = l + lmpml;
	    double cofl = sqrt(1.0 * (ltpo * (lmp + lmp + 1)));
	    summ = zpv[l30][ilmp30][0][0] /
	      (
	       dconjg(c1->rmi[l30][i40]) *
	       c1->rmi[lmp - 1][i40]
	       );
	    sume = zpv[l30][ilmp30][0][1] /
	      (
	       dconjg(c1->rmi[l30][i40]) *
	       c1->rei[lmp - 1][i40]
	       );
	    suem = zpv[l30][ilmp30][1][0] /
	      (
	       dconjg(c1->rei[l30][i40]) *
	       c1->rmi[lmp - 1][i40]
	       );
	    suee = zpv[l30][ilmp30][1][1] /
	      (
	       dconjg(c1->rei[l30][i40]) *
	       c1->rei[lmp - 1][i40]
	       );
	    sum += (cg1(lmpml, 0, l, -1) * (summ - sume - suem + suee) +
		    cg1(lmpml, 0, l, 1) * (summ + sume + suem + suee)) * cofl;
	  }
	}
      }
    }
    gaps[i40] = real(sum) * cofs;
  }
}

void cbf(int n, dcomplex z, int &nm, dcomplex *csj) {
  /*
   * FROM CSPHJY OF LIBRARY specfun
   *
   * ==========================================================
   *   Purpose: Compute spherical Bessel functions j
   *   Input :  z --- Complex argument of j
   *            n --- Order of j ( n = 0,1,2,... )
   *   Output:  csj(n+1) --- j
   *            nm --- Highest order computed
   *   Routines called:
   *            msta1 and msta2 for computing the starting
   *            point for backward recurrence
   * ==========================================================
   */
  double a0 = cabs(z);
  nm = n;
  if (a0 < 1.0e-60) {
    for (int k = 2; k <= n + 1; k++) {
      csj[k - 1] = 0.0;
    }
    csj[0] = 1.0 + 0.0 * I;
    return;
  }
  csj[0] = csin(z) / z;
  if (n == 0) {
    return;
  }
  csj[1] = (csj[0] - ccos(z)) / z;
  if (n == 1) {
    return;
  }
  dcomplex csa = csj[0];
  dcomplex csb = csj[1];
  int m = msta1(a0, 200);
  if (m < n) nm = m;
  else m = msta2(a0, n, 15);
  dcomplex cf0 = 0.0;
  dcomplex cf1 = 1.0e+00 - 100.0;
  dcomplex cf, cs;
  for (int k = m; k >= 0; k--) {
    cf = (2.0 * k + 3.0) * cf1 / z - cf0;
    if (k <= nm) csj[k] = cf;
    cf0 = cf1;
    cf1 = cf;
  }
  double abs_csa = cabs(csa);
  double abs_csb = cabs(csb);
  if (abs_csa > abs_csb) cs = csa / cf;
  else cs = csb / cf0;
  for (int k = 0; k <= nm; k++) {
    csj[k] = cs * csj[k];
  }
}

double cg1(int lmpml, int mu, int l, int m) {
  double result = 0.0;
  double xd, xn;
  if (lmpml == -1) { // Interpreted as GOTO 30
    xd = 2.0 * l * (2 * l - 1);
    if (mu == -1) {
      xn = 1.0 * (l - 1 - m) * (l - m);
    } else if (mu == 0) {
      xn = 2.0 * (l - m) * (l + m);
    } else if (mu == 1) {
      xn = 1.0 * (l - 1 + m) * (l + m);
    } else {
      throw 111; // Need an exception for unpredicted indices.
    }
    result = sqrt(xn / xd);
  } else if (lmpml == 0) { // Interpreted as GOTO 5
    bool goto10 = (m != 0) || (mu != 0);
    if (!goto10) {
      result = 0.0;
      return result;
    }
    if (mu != 0) {
      xd = 2.0 * l * (l + 1);
      if (mu == -1) {
	xn = 1.0 * (l - m) * (l + m + 1);
	result = -sqrt(xn / xd);
      } else if (mu == 1) { // mu > 0
	xn = 1.0 * (l + m) * (l - m + 1);
	result = sqrt(xn / xd);
      } else {
	throw 111; // Need an exception for unpredicted indices.
      }
    } else { // mu = 0
      xd = 1.0 * l * (l + 1);
      xn = -1.0 * m;
      result = xn / sqrt(xd);
    }
  } else if (lmpml == 1) { // Interpreted as GOTO 60
    xd = 2.0 * (l * 2 + 3) * (l + 1);
    if (mu == -1) {
      xn = 1.0 * (l + 1 + m) * (l + 2 + m);
      result = sqrt(xn / xd);
    } else if (mu == 0) {
      xn = 2.0 * (l + 1 - m) * (l + 1 + m);
      result = -sqrt(xn / xd);
    } else if (mu == 1) {
      xn = 1.0 * (l + 1 - m) * (l + 2 - m);
      result = sqrt(xn / xd);
    } else { // mu was not recognized.
      throw 111; // Need an exception for unpredicted indices.
    }
  } else { // lmpml was not recognized
    throw 111; // Need an exception for unpredicted indices.
  }
  return result;
}

void diel(int npntmo, int ns, int i, int ic, double vk, ParticleDescriptor *c1) {
  const double dif = c1->rc[i - 1][ns] - c1->rc[i - 1][ns - 1];
  const double half_step = 0.5 * dif / npntmo;
  double rr = c1->rc[i - 1][ns - 1];
  const dcomplex delta = c1->dc0[ic] - c1->dc0[ic - 1];
  const int kpnt = npntmo + npntmo;
  c1->ris[kpnt] = c1->dc0[ic];
  c1->dlri[kpnt] = 0.0 + 0.0 * I;
  const int i90 = i - 1;
  const int ns90 = ns - 1;
  const int ic90 = ic - 1;
  for (int np90 = 0; np90 < kpnt; np90++) {
    double ff = (rr - c1->rc[i90][ns90]) / dif;
    c1->ris[np90] = delta * ff * ff * (-2.0 * ff + 3.0) + c1->dc0[ic90];
    c1->dlri[np90] = 3.0 * delta * ff * (1.0 - ff) / (dif * vk * c1->ris[np90]);
    rr += half_step;
  }
}

void dme(
	 int li, int i, int npnt, int npntts, double vk, double exdc, double exri,
	 ParticleDescriptor *c1, int &jer, int &lcalc, dcomplex &arg, int last_conf
) {
  const int lipo = li + 1;
  const int lipt = li + 2;
  double *rfj = new double[lipt];
  double *rfn = new double[lipt];
  dcomplex cfj[lipt], fbi[lipt], fb[lipt], fn[lipt];
  dcomplex rmf[li], drmf[li], ref[li], dref[li];
  dcomplex dfbi, dfb, dfn, ccna, ccnb, ccnc, ccnd;
  dcomplex y1, dy1, y2, dy2, arin, cri;
  const dcomplex uim = 1.0 * I;
  int sph_index = (last_conf == 0) ? i : last_conf;
  jer = 0;
  int nstp = npnt - 1;
  int nstpts = npntts - 1;
  double sz = vk * c1->ros[sph_index - 1];
  c1->vsz[i - 1] = sz;
  double vkr1 = vk * c1->rc[sph_index - 1][0];
  int nsh = c1->nshl[sph_index - 1];
  c1->vkt[i - 1] = csqrt(c1->dc0[0]);
  arg = vkr1 * c1->vkt[i - 1];
  arin = arg;
  bool goto32 = false;
  if (imag(arg) != 0.0) {
    cbf(lipo, arg, lcalc, cfj);
    if (lcalc < lipo) {
      jer = 5;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    for (int j24 = 1; j24 <= lipt; j24++) fbi[j24 - 1] = cfj[j24 - 1];
    goto32 = true;
  }
  if (!goto32) {
    rbf(lipo, real(arg), lcalc, rfj);
    if (lcalc < lipo) {
      jer = 5;
      delete[] rfj;
      delete[] rfn;
      return;
    }
    for (int j30 = 1; j30 <= lipt; j30++) fbi[j30 - 1] = rfj[j30 - 1];
  }
  double arex = sz * exri;
  arg = arex;
  rbf(lipo, arex, lcalc, rfj);
  if (lcalc < lipo) {
    jer = 7;
    delete[] rfj;
    delete[] rfn;
    return;
  }
  rnf(lipo, arex, lcalc, rfn);
  if (lcalc < lipo) {
    jer = 8;
    delete[] rfj;
    delete[] rfn;
    return;
  }
  for (int j43 = 1; j43 <= lipt; j43++) {
    fb[j43 - 1] = rfj[j43 - 1];
    fn[j43 - 1] = rfn[j43 - 1];
  }
  if (nsh <= 1) {
    cri = c1->dc0[0] / exdc;
    for (int l60 = 1; l60 <= li; l60++) {
      int lpo = l60 + 1;
      int ltpo = lpo + l60;
      int lpt = lpo + 1;
      dfbi = ((1.0 * l60) * fbi[l60 - 1] - (1.0 * lpo) * fbi[lpt - 1]) * arin + fbi[lpo - 1] * (1.0 * ltpo);
      dfb = ((1.0 * l60) * fb[l60 - 1] - (1.0 * lpo) * fb[lpt - 1]) * arex + fb[lpo - 1] * (1.0 * ltpo);
      dfn = ((1.0 * l60) * fn[l60 - 1] - (1.0 * lpo) * fn[lpt - 1]) * arex + fn[lpo - 1] * (1.0 * ltpo);
      ccna = fbi[lpo - 1] * dfn;
      ccnb = fn[lpo - 1] * dfbi;
      ccnc = fbi[lpo - 1] * dfb;
      ccnd = fb[lpo - 1] * dfbi;
      c1->rmi[l60 - 1][i - 1] = 1.0 + uim * (ccna - ccnb) / (ccnc - ccnd);
      c1->rei[l60 - 1][i - 1] = 1.0 + uim * (cri * ccna - ccnb) / (cri * ccnc - ccnd);
    }
  } else { // nsh > 1
    int ic = 1;
    for (int l80 = 1; l80 <= li; l80++) {
      int lpo = l80 + 1;
      int ltpo = lpo + l80;
      int lpt = lpo + 1;
      int dltpo = ltpo;
      y1 = fbi[lpo - 1];
      dy1 = ((1.0 * l80) * fbi[l80 - 1] - (1.0 * lpo) * fbi[lpt - 1]) * c1->vkt[i - 1] / (1.0 * dltpo);
      y2 = y1;
      dy2 = dy1;
      ic = 1;
      for (int ns76 = 2; ns76 <= nsh; ns76++) {
	int nsmo = ns76 - 1;
	double vkr = vk * c1->rc[i - 1][nsmo - 1];
	if (ns76 % 2 != 0) {
	  ic += 1;
	  double step = 1.0 * nstp;
	  step = vk * (c1->rc[i - 1][ns76 - 1] - c1->rc[i - 1][nsmo - 1]) / step;
	  arg = c1->dc0[ic - 1];
	  rkc(nstp, step, arg, vkr, lpo, y1, y2, dy1, dy2);
	} else {
	  diel(nstpts, nsmo, i, ic, vk, c1);
	  double stepts = 1.0 * nstpts;
	  stepts = vk * (c1->rc[i - 1][ns76 - 1] - c1->rc[i - 1][nsmo - 1]) / stepts;
	  rkt(nstpts, stepts, vkr, lpo, y1, y2, dy1, dy2, c1);
	}
      }
      rmf[l80 - 1] = y1 * sz;
      drmf[l80 - 1] = dy1 * sz + y1;
      ref[l80 - 1] = y2 * sz;
      dref[l80 - 1] = dy2 * sz + y2;
    }
    cri = 1.0 + uim * 0.0;
    if (nsh % 2 != 0) cri = c1->dc0[ic - 1] / exdc;
    for (int l90 = 1; l90 <= li; l90++) {
      int lpo = l90 + 1;
      int ltpo = lpo + l90;
      int lpt = lpo + 1;
      dfb = ((1.0 * l90) * fb[l90 - 1] - (1.0 * lpo) * fb[lpt - 1]) * arex + fb[lpo - 1] * (1.0 * ltpo);
      dfn = ((1.0 * l90) * fn[l90 - 1] - (1.0 * lpo) * fn[lpt - 1]) * arex + fn[lpo - 1] * (1.0 * ltpo);
      ccna = rmf[l90 - 1] * dfn;
      ccnb = drmf[l90 - 1] * fn[lpo - 1] * (1.0 * sz * ltpo);
      ccnc = rmf[l90 - 1] * dfb;
      ccnd = drmf[l90 - 1] * fb[lpo -1] * (1.0 * sz * ltpo);
      c1->rmi[l90 - 1][i - 1] = 1.0 + uim *(ccna - ccnb) / (ccnc - ccnd);
      ccna = ref[l90 - 1] * dfn;
      ccnb = dref[l90 - 1] * fn[lpo - 1] * (1.0 * sz * ltpo);
      ccnc = ref[l90 - 1] * dfb;
      ccnd = dref[l90 - 1] * fb[lpo - 1] * (1.0 * sz * ltpo);
      c1->rei[l90 - 1][i - 1] = 1.0 + uim * (cri * ccna - ccnb) / (cri * ccnc - ccnd);
    }
  } // nsh <= 1 ?
  delete[] rfj;
  delete[] rfn;
  return;
}

double envj(int n, double x) {
  double result = 0.0;
  double xn;
  if (n == 0) {
    xn = 1.0e-100;
    result = 0.5 * log10(6.28 * xn) - xn * log10(1.36 * x / xn);
  } else {
    result = 0.5 * log10(6.28 * n) - n * log10(1.36 * x / n);
  }
  return result;
}

void mmulc(dcomplex *vint, double **cmullr, double **cmul) {
  double sm2 = real(vint[0]);
  double s24 = real(vint[1]);
  double d24 = imag(vint[1]);
  double sm4 = real(vint[5]);
  double s23 = real(vint[8]);
  double d32 = imag(vint[8]);
  double s34 = real(vint[9]);
  double d34 = imag(vint[9]);
  double sm3 = real(vint[10]);
  double s31 = real(vint[11]);
  double d31 = imag(vint[11]);
  double s21 = real(vint[12]);
  double d12 = imag(vint[12]);
  double s41 = real(vint[13]);
  double d14 = imag(vint[13]);
  double sm1 = real(vint[15]);
  cmullr[0][0] = sm2;
  cmullr[0][1] = sm3;
  cmullr[0][2] = -s23;
  cmullr[0][3] = -d32;
  cmullr[1][0] = sm4;
  cmullr[1][1] = sm1;
  cmullr[1][2] = -s41;
  cmullr[1][3] = -d14;
  cmullr[2][0] = -s24 * 2.0;
  cmullr[2][1] = -s31 * 2.0;
  cmullr[2][2] = s21 + s34;
  cmullr[2][3] = d34 + d12;
  cmullr[3][0] = -d24 * 2.0;
  cmullr[3][1] = -d31 * 2.0;
  cmullr[3][2] = d34 - d12;
  cmullr[3][3] = s21 - s34;
  cmul[0][0] = (sm2 + sm3 + sm4 + sm1) * 0.5;
  cmul[0][1] = (sm2 - sm3 + sm4 - sm1) * 0.5;
  cmul[0][2] = -s23 - s41;
  cmul[0][3] = -d32 - d14;
  cmul[1][0] = (sm2 + sm3 - sm4 - sm1) * 0.5;
  cmul[1][1] = (sm2 - sm3 - sm4 + sm1) * 0.5;
  cmul[1][2] = -s23 + s41;
  cmul[1][3] = -d32 + d14;
  cmul[2][0] = -s24 - s31;
  cmul[2][1] = -s24 + s31;
  cmul[2][2] = s21 + s34;
  cmul[2][3] = d34 + d12;
  cmul[3][0] = -d24 - d31;
  cmul[3][1] = -d24 + d31;
  cmul[3][2] = d34 - d12;
  cmul[3][3] = s21 - s34;
}

int msta1(double x, int mp) {
  int result = 0;
  double a0 = x;
  if (a0 < 0.0) a0 *= -1.0;
  int n0 = (int)(1.1 * a0) + 1;
  double f0 = envj(n0, a0) - mp;
  int n1 = n0 + 5;
  double f1 = envj(n1, a0) - mp;
  for (int it10 = 0; it10 < 20; it10++) {
    int nn = n1 - (int)((n1 - n0) / (1.0 - f0 / f1));
    double f = envj(nn, a0) - mp;
    int test_n = nn - n1;
    if (test_n < 0) test_n *= -1;
    if (test_n < 1) {
      return nn;
    }
    n0 = n1;
    f0 = f1;
    n1 = nn;
    f1 = f;
    result = nn;
  }
  return result;
}

int msta2(double x, int n, int mp) {
  int result = 0;
  double a0 = x;
  if (a0 < 0) a0 *= -1.0;
  double half_mp = 0.5 * mp;
  double ejn = envj(n, a0);
  double obj;
  int n0;
  if (ejn <= half_mp) {
    obj = 1.0 * mp;
    n0 = (int)(1.1 * a0) + 1;
  } else {
    obj = half_mp + ejn;
    n0 = n;
  }
  double f0 = envj(n0, a0) - obj;
  int n1 = n0 + 5;
  double f1 = envj(n1, a0) - obj;
  for (int it10 = 0; it10 < 20; it10 ++) {
    int nn = n1 - (int)((n1 - n0) / (1.0 - f0 / f1));
    double f = envj(nn, a0) - obj;
    int test_n = nn - n1;
    if (test_n < 0) test_n *= -1;
    if (test_n < 1) return (nn + 10);
    n0 = n1;
    f0 = f1;
    n1 = nn;
    f1 = f;
    result = nn + 10;
  }
  return result;
}

void orunve(double *u1, double *u2, double *u3, int iorth, double torth) {
  if (iorth <= 0) {
    double cp = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
    double abs_cp = cp;
    if (abs_cp < 0.0) abs_cp *= -1.0;
    if (iorth == 0 || abs_cp >= torth) {
      double fn = 1.0 / sqrt(1.0 - cp * cp);
      u3[0] = (u1[1] * u2[2] - u1[2] * u2[1]) * fn;
      u3[1] = (u1[2] * u2[0] - u1[0] * u2[2]) * fn;
      u3[2] = (u1[0] * u2[1] - u1[1] * u2[0]) * fn;
      return;
    }
  }
  u3[0] = u1[1] * u2[2] - u1[2] * u2[1];
  u3[1] = u1[2] * u2[0] - u1[0] * u2[2];
  u3[2] = u1[0] * u2[1] - u1[1] * u2[0];
}

void pwma(
	  double *up, double *un, dcomplex *ylm, int inpol, int lw,
	  int isq, ParticleDescriptor *c1
) {
  const double four_pi = 8.0 * acos(0.0);
  int is = isq;
  if (isq == -1) is = 0;
  int ispo = is + 1;
  int ispt = is + 2;
  int nlwm = lw * (lw + 2);
  int nlwmt = nlwm + nlwm;
  const double sqrtwi = 1.0 / sqrt(2.0);
  const dcomplex uim = 1.0 * I;
  dcomplex cm1 = 0.5 * (up[0] + up[1] * I);
  dcomplex cp1 = 0.5 * (up[0] - up[1] * I);
  double cz1 = up[2];
  dcomplex cm2 = 0.5 * (un[0] + un[1] * I);
  dcomplex cp2 = 0.5 * (un[0] - un[1] * I);
  double cz2 = un[2];
  for (int l20 = 0; l20 < lw; l20++) {
    int l = l20 + 1;
    int lf = l + 1;
    int lftl = lf * l;
    double x = 1.0 * lftl;
    dcomplex cl = four_pi / sqrt(x) + 0.0 * I;
    for (int powi = 1; powi <= l; powi++) cl *= uim;
    int mv = l + lf;
    int m = -lf;
    for (int mf20 = 0; mf20 < mv; mf20++) {
      m += 1;
      int k = lftl + m;
      x = 1.0 * (lftl - m * (m + 1));
      double cp = sqrt(x);
      x = 1.0 * (lftl - m * (m - 1));
      double cm = sqrt(x);
      double cz = 1.0 * m;
      c1->w[k - 1][ispo - 1] = dconjg(
				      cp1 * cp * ylm[k + 1] +
				      cm1 * cm * ylm[k - 1] +
				      cz1 * cz * ylm[k]
				      ) * cl;
      c1->w[k - 1][ispt - 1] = dconjg(
				      cp2 * cp * ylm[k + 1] +
				      cm2 * cm * ylm[k - 1] +
				      cz2 * cz * ylm[k]
				      ) * cl;
    }
  }
  for (int k30 = 0; k30 < nlwm; k30++) {
    int i = k30 + nlwm;
    c1->w[i][ispo - 1] = uim * c1->w[k30][ispt - 1];
    c1->w[i][ispt - 1] = -uim * c1->w[k30][ispo - 1];
  }
  if (inpol != 0) {
    for (int k40 = 0; k40 < nlwm; k40++) {
      int i = k40 + nlwm;
      dcomplex cc1 = sqrtwi * (c1->w[k40][ispo - 1] + uim * c1->w[k40][ispt - 1]);
      dcomplex cc2 = sqrtwi * (c1->w[k40][ispo - 1] - uim * c1->w[k40][ispt - 1]);
      c1->w[k40][ispo - 1] = cc2;
      c1->w[i][ispo - 1] = -cc2;
      c1->w[k40][ispt - 1] = cc1;
      c1->w[i][ispt - 1] = cc1;
    }
  } else {
    if (isq == 0) {
      return;
    }
  }
  if (isq != 0) {
    for (int i50 = 0; i50 < 2; i50++) {
      int ipt = i50 + 2;
      int ipis = i50 + is;
      for (int k50 = 0; k50 < nlwmt; k50++) {
	c1->w[k50][ipt] = dconjg(c1->w[k50][ipis]);
      }
    }
  }
}

void rabas(
	   int inpol, int li, int nsph, ParticleDescriptor *c1, double **tqse, dcomplex **tqspe,
	   double **tqss, dcomplex **tqsps
) {
  dcomplex cc0 = 0.0 + 0.0 * I;
  dcomplex uim = 0.0 + 1.0 * I;
  double two_pi = 4.0 * acos(0.0);
  for (int i80 = 0; i80 < nsph; i80++) {
    int i = i80 + 1;
    if(c1->iog[i80] >= i) {
      tqse[0][i80] = 0.0;
      tqse[1][i80] = 0.0;
      tqspe[0][i80] = cc0;
      tqspe[1][i80] = cc0;
      tqss[0][i80] = 0.0;
      tqss[1][i80] = 0.0;
      tqsps[0][i80] = cc0;
      tqsps[1][i80] = cc0;
      for (int l70 = 0; l70 < li; l70++) {
	int l = l70 + 1;
	double fl = 1.0 * (l + l + 1);
	dcomplex rm = 1.0 / c1->rmi[l70][i80];
	double rmm = real(rm * dconjg(rm));
	dcomplex re = 1.0 / c1->rei[l70][i80];
	double rem = real(re * dconjg(re));
	if (inpol == 0) {
	  dcomplex pce = ((rm + re) * uim) * fl;
	  dcomplex pcs = ((rmm + rem) * fl) * uim;
	  tqspe[0][i80] -= pce;
	  tqspe[1][i80] += pce;
	  tqsps[0][i80] -= pcs;
	  tqsps[1][i80] += pcs;
	} else {
	  double ce = real(rm + re) * fl;
	  double cs = (rmm + rem) * fl;
	  tqse[0][i80] -= ce;
	  tqse[1][i80] += ce;
	  tqss[0][i80] -= cs;
	  tqss[1][i80] += cs;
	}
      }
      if (inpol == 0) {
	tqspe[0][i80] *= two_pi;
	tqspe[1][i80] *= two_pi;
	tqsps[0][i80] *= two_pi;
	tqsps[1][i80] *= two_pi;
      } else {
	tqse[0][i80] *= two_pi;
	tqse[1][i80] *= two_pi;
	tqss[0][i80] *= two_pi;
	tqss[1][i80] *= two_pi;
      }
    }
  }
}

void rbf(int n, double x, int &nm, double sj[]) {
  /*
   * FROM SPHJ OF LIBRARY specfun
   *
   * ==========================================================
   *   Purpose: Compute spherical Bessel functions j
   *   Input :  x --- Argument of j
   *            n --- Order of j ( n = 0,1,2,... )
   *   Output:  sj(n+1) --- j
   *            nm --- Highest order computed
   *   Routines called:
   *            msta1 and msta2 for computing the starting
   *            point for backward recurrence
   * ==========================================================
   */
  double a0 = x;
  if (a0 < 0.0) a0 *= -1.0;
  nm = n;
  if (a0 < 1.0e-60) {
    for (int k = 1; k <= n; k++)
      sj[k] = 0.0;
    sj[0] = 1.0;
    return;
  }
  sj[0] = sin(x) / x;
  if (n == 0) {
    return;
  }
  sj[1] = (sj[0] - cos(x)) / x;
  if (n == 1) {
    return;
  }
  double sa = sj[0];
  double sb = sj[1];
  int m = msta1(a0, 200);
  if (m < n) nm = m;
  else m = msta2(a0, n, 15);
  double f0 = 0.0;
  double f1 = 1.0e+00 - 100.0;
  double f;
  for (int k = m; k >= 0; k--) {
    f = (2.0 * k +3.0) * f1 / x - f0;
    if (k <= nm) sj[k] = f;
    f0 = f1;
    f1 = f;
  }
  double cs;
  double abs_sa = (sa < 0.0) ? -sa : sa;
  double abs_sb = (sb < 0.0) ? -sb : sb;
  if (abs_sa > abs_sb) cs = sa / f;
  else cs = sb / f0;
  for (int k = 0; k <= nm ; k++) {
    sj[k] = cs * sj[k];
  }
}

void rkc(
	 int npntmo, double step, dcomplex dcc, double &x, int lpo,
	 dcomplex &y1, dcomplex &y2, dcomplex &dy1, dcomplex &dy2
) {
  dcomplex cy1, cdy1, c11, cy23, yc2, c12, c13;
  dcomplex cy4, yy, c14, c21, c22, c23, c24;
  double half_step = 0.5 * step;
  double cl = 1.0 * lpo * (lpo - 1);
  for (int ipnt60 = 0; ipnt60 < npntmo; ipnt60++) {
    cy1 = cl / (x * x) - dcc;
    cdy1 = -2.0 / x;
    c11 = (cy1 * y1 + cdy1 * dy1) * step;
    double xh = x + half_step;
    cy23 = cl / (xh * xh) - dcc;
    double cdy23 = -2.0 / xh;
    yc2 = y1 + dy1 * half_step;
    c12 = (cy23 * yc2 + cdy23 * (dy1 + 0.5 * c11)) * step;
    c13 = (cy23 * (yc2 + 0.25 * c11 * step) + cdy23 * (dy1 + 0.5 * c12)) * step;
    double xn = x + step;
    cy4 = cl / (xn * xn) - dcc;
    double cdy4 = -2.0 / xn;
    yy = y1 + dy1 * step;
    c14 = (cy4 * (yy + 0.5 * c12 * step) + cdy4 * (dy1 + c13)) * step;
    y1 = yy + (c11 + c12 + c13) * step / 6.0;
    dy1 += (0.5 * c11 + c12 + c13 + 0.5 * c14) / 3.0;
    c21 = (cy1 * y2 + cdy1 * dy2) * step;
    yc2 = y2 + dy2 * half_step;
    c22 = (cy23 * yc2 + cdy23 * (dy2 + 0.5 * c21)) * step;
    c23 = (cy23 * (yc2 + 0.25 * c21 * step) + cdy23 * (dy2 + 0.5 * c22)) * step;
    yy = y2 + dy2 * step;
    c24 = (cy4 * (yc2 + 0.5 * c22 * step) + cdy4 * (dy2 + c23)) * step;
    y2 = yy + (c21 + c22 + c23) * step / 6.0;
    dy2 += (0.5 * c21 + c22 + c23 + 0.5 * c24) / 3.0;
    x = xn;
  }
}

void rkt(
	 int npntmo, double step, double &x, int lpo, dcomplex &y1,
	 dcomplex &y2, dcomplex &dy1, dcomplex &dy2, ParticleDescriptor *c1
) {
  dcomplex cy1, cdy1, c11, cy23, cdy23, yc2, c12, c13;
  dcomplex cy4, cdy4, yy, c14, c21, c22, c23, c24;
  double half_step = 0.5 * step;
  double cl = 1.0 * lpo * (lpo - 1);
  for (int ipnt60 = 0; ipnt60 < npntmo; ipnt60++) {
    int ipnt = ipnt60 + 1;
    int jpnt = ipnt + ipnt - 1;
    int jpnt60 = jpnt - 1;
    cy1 = cl / (x * x) - c1->ris[jpnt60];
    cdy1 = -2.0 / x;
    c11 = (cy1 * y1 + cdy1 * dy1) * step;
    double xh = x + half_step;
    int jpntpo = jpnt + 1;
    cy23 = cl / (xh * xh) - c1->ris[jpnt];
    cdy23 = -2.0 / xh;
    yc2 = y1 + dy1 * half_step;
    c12 = (cy23 * yc2 + cdy23 * (dy1 + 0.5 * c11)) * step;
    c13= (cy23 * (yc2 + 0.25 * c11 *step) + cdy23 * (dy1 + 0.5 * c12)) * step;
    double xn = x + step;
    //int jpntpt = jpnt + 2;
    cy4 = cl / (xn * xn) - c1->ris[jpntpo];
    cdy4 = -2.0 / xn;
    yy = y1 + dy1 * step;
    c14 = (cy4 * (yy + 0.5 * c12 * step) + cdy4 * (dy1 + c13)) * step;
    y1= yy + (c11 + c12 + c13) * step / 6.0;
    dy1 += (0.5 * c11 + c12 + c13 + 0.5 * c14) /3.0;
    cy1 -= cdy1 * c1->dlri[jpnt60];
    cdy1 += 2.0 * c1->dlri[jpnt60];
    c21 = (cy1 * y2 + cdy1 * dy2) * step;
    cy23 -= cdy23 * c1->dlri[jpnt];
    cdy23 += 2.0 * c1->dlri[jpnt];
    yc2 = y2 + dy2 * half_step;
    c22 = (cy23 * yc2 + cdy23 * (dy2 + 0.5 * c21)) * step;
    c23 = (cy23 * (yc2 + 0.25 * c21 * step) + cdy23 * (dy2 + 0.5 * c22)) * step;
    cy4 -= cdy4 * c1->dlri[jpntpo];
    cdy4 += 2.0 * c1->dlri[jpntpo];
    yy = y2 + dy2 * step;
    c24 = (cy4 * (yc2 + 0.5 * c22 * step) + cdy4 * (dy2 + c23)) * step;
    y2 = yy + (c21 + c22 + c23) * step / 6.0;
    dy2 += (0.5 * c21 + c22 + c23 + 0.5 * c24) / 3.0;
    x = xn;
  }
}

void rnf(int n, double x, int &nm, double *sy) {
  /*
   * FROM SPHJY OF LIBRARY specfun
   *
   * ==========================================================
   *   Purpose: Compute spherical Bessel functions y
   *   Input :  x --- Argument of y ( x > 0 )
   *            n --- Order of y ( n = 0,1,2,... )
   *   Output:  sy(n+1) --- y
   *            nm --- Highest order computed
   * ==========================================================
   */
  if (x < 1.0e-60) {
    for (int k = 0; k <= n; k++)
      sy[k] = -1.0e300;
    return;
  }
  sy[0] = -1.0 * cos(x) / x;
  if (n == 0) {
    return;
  }
  sy[1] = (sy[0] - sin(x)) / x;
  if (n == 1) {
    return;
  }
  double f0 = sy[0];
  double f1 = sy[1];
  double f;
  for (int k = 2; k <= n; k++) {
    f = (2.0 * k - 1.0) * f1 / x - f0;
    sy[k] = f;
    double abs_f = f;
    if (abs_f < 0.0) abs_f *= -1.0;
    if (abs_f >= 1.0e300) {
      nm = k;
      break;
    }
    f0 = f1;
    f1 = f;
    nm = k;
  }
  return;
}

void sphar(
	   double cosrth, double sinrth, double cosrph, double sinrph,
	   int ll, dcomplex *ylm
) {
  const int rmp_size = ll;
  const int plegn_size = (ll + 1) * ll / 2 + ll + 1;
  double sinrmp[rmp_size], cosrmp[rmp_size], plegn[plegn_size];
  double four_pi = 8.0 * acos(0.0);
  double pi4irs = 1.0 / sqrt(four_pi);
  double x = cosrth;
  double y = sinrth;
  if (y < 0.0) y *= -1.0;
  double cllmo = 3.0;
  double cll = 1.5;
  double ytol = y;
  plegn[0] = 1.0;
  plegn[1] = x * sqrt(cllmo);
  plegn[2] = ytol * sqrt(cll);
  sinrmp[0] = sinrph;
  cosrmp[0] = cosrph;
  if (ll >= 2) {
    int k = 3;
    for (int l20 = 2; l20 <= ll; l20++) {
      int lmo = l20 - 1;
      int ltpo = l20 + l20 + 1;
      int ltmo = ltpo - 2;
      int lts = ltpo * ltmo;
      double cn = 1.0 * lts;
      for (int mpo10 = 1; mpo10 <= lmo; mpo10++) {
	int m = mpo10 - 1;
	int mpopk = mpo10 + k;
	int ls = (l20 + m) * (l20 - m);
	double cd = 1.0 * ls;
	double cnm = 1.0 * ltpo * (lmo + m) * (l20 - mpo10);
	double cdm = 1.0 * ls * (ltmo - 2);
	plegn[mpopk - 1] = plegn[mpopk - l20 - 1] * x * sqrt(cn / cd) -
	  plegn[mpopk - ltmo - 1] * sqrt(cnm / cdm);
      }
      int lpk = l20 + k;
      double cltpo = 1.0 * ltpo;
      plegn[lpk - 1] = plegn[k - 1] * x * sqrt(cltpo);
      k = lpk + 1;
      double clt = 1.0 * (ltpo - 1);
      cll *= (cltpo / clt);
      ytol *= y;
      plegn[k - 1] = ytol * sqrt(cll);
      sinrmp[l20 - 1] = sinrph * cosrmp[lmo - 1] + cosrph * sinrmp[lmo - 1];
      cosrmp[l20 - 1] = cosrph * cosrmp[lmo - 1] - sinrph * sinrmp[lmo - 1];
    } // end l20 loop
  }
  // label 30
  int l = 0;
  int m, k, l0y, l0p, lmy, lmp;
  double save;
 label40:
  m = 0;
  k = l * (l + 1);
  l0y = k + 1;
  l0p = k / 2 + 1;
  ylm[l0y - 1] = pi4irs * plegn[l0p - 1];
  goto label45;
 label44:
  lmp = l0p + m;
  save = pi4irs * plegn[lmp - 1];
  lmy = l0y + m;
  ylm[lmy - 1] = save * (cosrmp[m - 1] + sinrmp[m - 1] * I);
  if (m % 2 != 0) ylm[lmy - 1] *= -1.0;
  lmy = l0y - m;
  ylm[lmy - 1] = save * (cosrmp[m - 1] - sinrmp[m - 1] * I);
 label45:
  if (m >= l) goto label47;
  m += 1;
  goto label44;
 label47:
  if (l >= ll) return;
  l += 1;
  goto label40;
}

void sscr0(dcomplex &tfsas, int nsph, int lm, double vk, double exri, ParticleDescriptor *c1) {
  dcomplex sum21, rm, re, csam;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  const double exdc = exri * exri;
  double ccs = 4.0 * acos(0.0) / (vk * vk);
  double cccs = ccs / exdc;
  csam = -(ccs / (exri * vk)) * (0.0 + 0.5 * I);
  tfsas = cc0;
  for (int i12 = 0; i12 < nsph; i12++) {
    int i = i12 + 1;
    int iogi = c1->iog[i12];
    if (iogi >= i) {
      double sums = 0.0;
      dcomplex sum21 = cc0;
      for (int l10 = 0; l10 < lm; l10++) {
	int l = l10 + 1;
	double fl = 1.0 + l + l;
	rm = 1.0 / c1->rmi[l10][i12];
	re = 1.0 / c1->rei[l10][i12];
	dcomplex rm_cnjg = dconjg(rm);
	dcomplex re_cnjg = dconjg(re);
	sums += real(rm_cnjg * rm + re_cnjg * re) * fl;
	sum21 += (rm + re) * fl;
      }
      sum21 *= -1.0;
      double scasec = cccs * sums;
      double extsec = -cccs * real(sum21);
      double abssec = extsec - scasec;
      c1->sscs[i12] = scasec;
      c1->sexs[i12] = extsec;
      c1->sabs[i12] = abssec;
      double gcss = c1->gcsv[i12];
      c1->sqscs[i12] = scasec / gcss;
      c1->sqexs[i12] = extsec / gcss;
      c1->sqabs[i12] = abssec / gcss;
      c1->fsas[i12] = sum21 * csam;
    }
    tfsas += c1->fsas[iogi - 1];
  }
}

void sscr2(int nsph, int lm, double vk, double exri, ParticleDescriptor *c1) {
  dcomplex s11, s21, s12, s22, rm, re, csam, cc;
  const dcomplex cc0 = 0.0 + 0.0 * I;
  double ccs = 1.0 / (vk * vk);
  csam = -(ccs / (exri * vk)) * (0.0 + 0.5 * I);
  const double pigfsq = 64.0 * acos(0.0) * acos(0.0);
  double cfsq = 4.0 / (pigfsq * ccs * ccs);
  int nlmm = lm * (lm + 2);
  for (int i14 = 0; i14 < nsph; i14++) {
    int i = i14 + 1;
    int iogi = c1->iog[i14];
    if (iogi >= i) {
      int k = 0;
      s11 = cc0;
      s21 = cc0;
      s12 = cc0;
      s22 = cc0;
      for (int l10 = 0; l10 < lm; l10++) {
	int l = l10 + 1;
	rm = 1.0 / c1->rmi[l10][i14];
	re = 1.0 / c1->rei[l10][i14];
	int ltpo = l + l + 1;
	for (int im10 = 0; im10 < ltpo; im10++) {
	  k += 1;
	  int ke = k + nlmm;
	  s11 = s11 - c1->w[k - 1][2] * c1->w[k - 1][0] * rm - c1->w[ke - 1][2] * c1->w[ke - 1][0] * re;
	  s21 = s21 - c1->w[k - 1][3] * c1->w[k - 1][0] * rm - c1->w[ke - 1][3] * c1->w[ke - 1][0] * re;
	  s12 = s12 - c1->w[k - 1][2] * c1->w[k - 1][1] * rm - c1->w[ke - 1][2] * c1->w[ke - 1][1] * re;
	  s22 = s22 - c1->w[k - 1][3] * c1->w[k - 1][1] * rm - c1->w[ke - 1][3] * c1->w[ke - 1][1] * re;
	}
      }
      c1->sas[i14][0][0] = s11 * csam;
      c1->sas[i14][1][0] = s21 * csam;
      c1->sas[i14][0][1] = s12 * csam;
      c1->sas[i14][1][1] = s22 * csam;
    }
  } // loop i14
  for (int i24 = 0; i24 < nsph; i24++) {
    int i = i24 + 1;
    int iogi = c1->iog[i24];
    if (iogi >= i) {
      int j = 0;
      for (int ipo1 = 0; ipo1 < 2; ipo1++) {
	for (int jpo1 = 0; jpo1 < 2; jpo1++) {
	  dcomplex cc = dconjg(c1->sas[i24][jpo1][ipo1]);
	  for (int ipo2 = 0; ipo2 < 2; ipo2++) {
	    for (int jpo2 = 0; jpo2 < 2; jpo2++) {
	      c1->vints[i24][j++] = c1->sas[i24][jpo2][ipo2] * cc * cfsq;
	    }
	  }
	}
      }
    }
  }
}

void thdps(int lm, double ****zpv) {
  for (int l15 = 0; l15 < lm; l15++) {
    int l = l15 + 1;
    double xd = 1.0 * l * (l + 1);
    double zp = -1.0 / sqrt(xd);
    zpv[l15][1][0][1] = zp;
    zpv[l15][1][1][0] = zp;
  }
  if (lm != 1) {
    for (int l20 = 1; l20 < lm; l20++) {
      int l = l20 + 1;
      double xn = 1.0 * (l - 1) * (l + 1);
      double xd = 1.0 * l * (l + l + 1);
      double zp = sqrt(xn / xd);
      zpv[l20][0][0][0] = zp;
      zpv[l20][0][1][1] = zp;
    }
    int lmmo = lm - 1;
    for (int l25 = 0; l25 < lmmo; l25++) {
      int l = l25 + 1;
      double xn = 1.0 * l * (l + 2);
      double xd = (l + 1) * (l + l + 1);
      double zp = -1.0 * sqrt(xn / xd);
      zpv[l25][2][0][0] = zp;
      zpv[l25][2][1][1] = zp;
    }
  }
}

void upvmp(
	   double thd, double phd, int icspnv, double &cost, double &sint,
	   double &cosp, double &sinp, double *u, double *up, double *un
) {
  double half_pi = acos(0.0);
  double rdr = half_pi / 90.0;
  double th = thd * rdr;
  double ph = phd * rdr;
  cost = cos(th);
  sint = sin(th);
  cosp = cos(ph);
  sinp = sin(ph);
  u[0] = cosp * sint;
  u[1] = sinp * sint;
  u[2] = cost;
  up[0] = cosp * cost;
  up[1] = sinp * cost;
  up[2] = -sint;
  un[0] = -sinp;
  un[1] = cosp;
  un[2] = 0.0;
  if (icspnv != 0) {
    up[0] *= -1.0;
    up[1] *= -1.0;
    up[2] *= -1.0;
    un[0] *= -1.0;
    un[1] *= -1.0;
  }
}

void upvsp(
	   double *u, double *upmp, double *unmp, double *us, double *upsmp, double *unsmp,
	   double *up, double *un, double *ups, double *uns, double *duk, int &isq,
	   int &ibf, double &scand, double &cfmp, double &sfmp, double &cfsp, double &sfsp
) {
  double rdr = acos(0.0) / 90.0;
  double small = 1.0e-6;
  isq = 0;
  scand = u[0] * us[0] + u[1] * us[1] + u[2] * us[2];
  double abs_scand = (scand >= 1.0) ? scand - 1.0 : 1.0 - scand;
  if (abs_scand >= small) {
    abs_scand = scand + 1.0;
    if (abs_scand < 0.0) abs_scand *= -1.0;
    if (abs_scand >= small) {
      scand = acos(scand) / rdr;
      duk[0] = u[0] - us[0];
      duk[1] = u[1] - us[1];
      duk[2] = u[2] - us[2];
      ibf = 0;
    } else { // label 15
      scand = 180.0;
      duk[0] = 2.0 * u[0];
      duk[1] = 2.0 * u[1];
      duk[2] = 2.0 * u[2];
      ibf = 1;
      ups[0] = -upsmp[0];
      ups[1] = -upsmp[1];
      ups[2] = -upsmp[2];
      uns[0] = -unsmp[0];
      uns[1] = -unsmp[1];
      uns[2] = -unsmp[2];
    }
  } else { // label 10
    scand = 0.0;
    duk[0] = 0.0;
    duk[1] = 0.0;
    duk[2] = 0.0;
    ibf = -1;
    isq = -1;
    ups[0] = upsmp[0];
    ups[1] = upsmp[1];
    ups[2] = upsmp[2];
    uns[0] = unsmp[0];
    uns[1] = unsmp[1];
    uns[2] = unsmp[2];
  }
  if (ibf == -1 || ibf == 1) { // label 20
    up[0] = upmp[0];
    up[1] = upmp[1];
    up[2] = upmp[2];
    un[0] = unmp[0];
    un[1] = unmp[1];
    un[2] = unmp[2];
  } else { // label 25
    orunve(u, us, un, -1, small);
    uns[0] = un[0];
    uns[1] = un[1];
    uns[2] = un[2];
    orunve(un, u, up, 1, small);
    orunve(uns, us, ups, 1, small);
  }
  // label 85
  cfmp = upmp[0] * up[0] + upmp[1] * up[1] + upmp[2] * up[2];
  sfmp = unmp[0] * up[0] + unmp[1] * up[1] + unmp[2] * up[2];
  cfsp = ups[0] * upsmp[0] + ups[1] * upsmp[1] + ups[2] * upsmp[2];
  sfsp = uns[0] * upsmp[0] + uns[1] * upsmp[1] + uns[2] * upsmp[2];
}

void wmamp(
	   int iis, double cost, double sint, double cosp, double sinp, int inpol,
	   int lm, int idot, int nsph, double *arg, double *u, double *up,
	   double *un, ParticleDescriptor *c1
) {
  const int ylm_size = (lm + 1) * (lm + 1) + 1;
  dcomplex *ylm = new dcomplex[ylm_size];
  const int nlmp = lm * (lm + 2) + 2;
  ylm[nlmp - 1] = 0.0 + 0.0 * I;
  if (idot != 0) {
    if (idot != 1) {
      for (int n40 = 0; n40 < nsph; n40++) {
	arg[n40] = u[0] * c1->rxx[n40] + u[1] * c1->ryy[n40] + u[2] * c1->rzz[n40];
      }
    } else {
      for (int n50 = 0; n50 < nsph; n50++) {
	arg[n50] = c1->rzz[n50];
      }
    }
    if (iis == 2) {
      for (int n60 = 0; n60 < nsph; n60++) arg[n60] *= -1;
    }
  }
  sphar(cost, sint, cosp, sinp, lm, ylm);
  pwma(up, un, ylm, inpol, lm, iis, c1);
  delete[] ylm;
}

void wmasp(
	   double cost, double sint, double cosp, double sinp, double costs, double sints,
	   double cosps, double sinps, double *u, double *up, double *un, double *us,
	   double *ups, double *uns, int isq, int ibf, int inpol, int lm, int idot,
	   int nsph, double *argi, double *args, ParticleDescriptor *c1
) {
  const int ylm_size = (lm + 1) * (lm + 1) + 1;
  dcomplex *ylm = new dcomplex[ylm_size];
  const int nlmp = lm * (lm + 2) + 2;
  ylm[nlmp - 1] = 0.0 + 0.0 * I;
  if (idot != 0) {
    if (idot != 1) {
      for (int n40 = 0; n40 < nsph; n40++) {
	argi[n40] = u[0] * c1->rxx[n40] + u[1] * c1->ryy[n40] + u[2] * c1->rzz[n40];
	if (ibf != 0) {
	  args[n40] = argi[n40] * ibf;
	} else {
	  args[n40] = -1.0 * (us[0] * c1->rxx[n40] + us[1] * c1->ryy[n40] + us[2] * c1->rzz[n40]);
	}
      }
    } else { // label 50
      for (int n60 = 0; n60 < nsph; n60++) {
	argi[n60] = cost * c1->rzz[n60];
	if (ibf != 0) {
	  args[n60] = argi[n60] * ibf;
	} else {
	  args[n60] = -costs * c1->rzz[n60];
	}
      }
    }
  }
  sphar(cost, sint, cosp, sinp, lm, ylm);
  pwma(up, un, ylm, inpol, lm, isq, c1);
  if (ibf >= 0) {
    sphar(costs, sints, cosps, sinps, lm, ylm);
    pwma(ups, uns, ylm, inpol, lm, 2, c1);
  }
  delete[] ylm;
}
