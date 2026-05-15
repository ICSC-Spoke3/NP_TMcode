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

/*! \file magma_calls.cpp
 *
 * \brief Implementation of the interface with MAGMA libraries.
 */

#ifdef USE_MAGMA

#include <cstdio>
#include <cstring>

#ifdef USE_TARGET_OFFLOAD
#include "omp.h"
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
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

#ifndef INCLUDE_MAGMA_CALLS_H_
#include "../include/magma_calls.h"
#endif

#ifdef DEBUG_REFINE
#include <hdf5.h>
#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif // INCLUDE_ERRORS_H_
#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif // INCLUDE_LIST_H_
#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif // INCLUDE_FILE_IO_H_
#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif // INCLUDE_UTILS_H_

#endif // DEBUG_REFINE

using namespace std;

#define FOUR_PI 12.566370614359172

#ifdef USE_TARGET_OFFLOAD

#pragma omp begin declare target device_type(any)
inline magmaDoubleComplex magma_zadd(const magmaDoubleComplex& a, const magmaDoubleComplex& b) {
  return MAGMA_Z_MAKE(a.x+b.x,a.y+b.y);
}
#pragma omp end declare target

#pragma omp begin declare target device_type(any)
inline magmaDoubleComplex magma_zprod(const magmaDoubleComplex& a, const magmaDoubleComplex& b) {
  return MAGMA_Z_MAKE(
    a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x
  );
}
#pragma omp end declare target

#pragma omp begin declare target device_type(any)
inline magmaDoubleComplex magma_zscale(const magmaDoubleComplex a, const double b) {
  return MAGMA_Z_MAKE(a.x*b,a.y*b);
}
#pragma omp end declare target

#pragma omp begin declare target device_type(any)
double magma_cgev(int ipamo, int mu, int l, int m) {
  double xd, xn;
  if (ipamo == 0) {
    if (mu != 0) {
      xd = 2.0 * l * (l + 1);
      xn = (mu < 0) ? 1.0 * (l + m) * (l - m + 1)
	: 1.0 * (l - m) * (l + m + 1);
      double result = sqrt(xn / xd);
      return (mu <= 0) ? result : -result;
    } else { // mu == 0
      return -1.0 * m / sqrt(1.0 * (l + 1) * l);
    }
    return 0.0;
  } else { // ipamo != 0
    xd = 2.0 * l * (l * 2 - 1);
    xn = (mu == 0) ? 2.0 * (l - m) * (l + m)
      : ((mu < 0) ? 1.0 * (l - 1 + m) * (l + m) : 1.0 * (l - 1 - m) * (l - m));
    return sqrt(xn / xd);
  }
}
#pragma omp end declare target

magma_int_t magma_cms(
  magmaDoubleComplex *vec_am, ParticleDescriptor *c1, double *rxx, double *ryy, double *rzz,
  int *ind3j, double *v3j0, magmaDoubleComplex *vh, magmaDoubleComplex *vyhj,
  magmaDoubleComplex *vj0, magmaDoubleComplex *vyj0, magmaDoubleComplex *rmi,
  magmaDoubleComplex* rei, int device_id
) {
  const magmaDoubleComplex cz0 = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex cz1 = MAGMA_Z_MAKE(1.0, 0.0); 
  const int nsph = c1->nsph;
  const int nsphmo = nsph - 1;
  const int ncou = nsph * nsphmo;
  const int li = c1->li;
  const int litpo = li + li + 1;
  const int litpos = litpo * litpo;
  const int le = c1->le;
  const int lmtpo = c1->lmtpo;
  const int lmtpos = c1->lmtpos;
  const int rsize = li * nsph;
  const magma_int_t max_litpo = 2 * li + 1;
  const magma_int_t nlim = li * (li + 2);
  const magma_int_t ndi = nsph * nlim;
  const magma_int_t ndit = 2 * ndi;
  const magma_int_t size = ndit * ndit;
  const magma_int_t num_pairs = (nsph * (nsph - 1)) / 2;
  const magma_int_t total_iters = num_pairs * li * max_litpo * li * max_litpo;
  magmaDoubleComplex vj = MAGMA_Z_MAKE(real(c1->vj), imag(c1->vj));

  int lut_n1[num_pairs];
  int lut_n2[num_pairs];
  for (int n1 = 1; n1 <= nsphmo; n1++) {
    for (int n2 = n1 + 1; n2 <= nsph; n2++) {
      int nbl = ((n1 - 1) * (2 * nsph - n1)) / 2 + n2 - n1;
      lut_n1[nbl - 1] = n1;
      lut_n2[nbl - 1] = n2;
    }
  }

#pragma omp target teams distribute parallel for \
  map(to: lut_n1[0:num_pairs], lut_n2[0:num_pairs]) \
  device(device_id)
  for (magma_int_t iter = 0; iter < total_iters; ++iter) {
    magma_int_t t = iter;
    double rac3j[MAX_GPU_STACK_DOUBLE]; // MAX_GPU_STACK_DOUBLE = 128

    // Index calculation (from innermost to outermost)
    int im2 = (t % max_litpo) + 1;
    t /= max_litpo;
    int l2  = (t % li) + 1;
    t /= li;
    int im1 = (t % max_litpo) + 1;
    t /= max_litpo;
    int l1  = (t % li) + 1;
    t /= li;
    int nbl_idx = (t % num_pairs); // 0-based for look-up tables
    // Look-up values
    int n1 = lut_n1[nbl_idx];
    int n2 = lut_n2[nbl_idx];
    int nbl = nbl_idx + 1;
    // Ranges of magnetic quantum numbers
    int l1po = l1 + 1;
    int l2po = l2 + 1;
    int l1tpo = 2 * l1 + 1;
    int l2tpo = 2 * l2 + 1;

    int il1 = l1po * l1;
    int in1 = (n1 - 1) * nlim;
    int in2 = (n2 - 1) * nlim;
    int il2 = l2po * l2;
    magmaDoubleComplex rsh = ((l2 + l1) % 2 == 0) ? MAGMA_Z_MAKE(1.0, 0.0) : MAGMA_Z_MAKE(-1.0, 0.0);
    magmaDoubleComplex rsk = MAGMA_Z_NEGATE(rsh);
    bool is_valid_iter = (im1 <= l1tpo && im2 <= l2tpo);
    // Index 1 magnetic quantum numbers
    int m1 = (is_valid_iter) ? im1 - l1po : 0;
    int ilm1 = (is_valid_iter) ? il1 + m1 : 0;
    int ilm1e = (is_valid_iter) ? ilm1 + ndi : 0;
    int i1 = (is_valid_iter) ? in1 + ilm1 : 0;
    int i1e = (is_valid_iter) ? in1 + ilm1e : 0;
    int j1 = (is_valid_iter) ? in2 + ilm1 : 0;
    int j1e = (is_valid_iter) ? in2 + ilm1e : 0;
    // End of index 1 magnetic quantum numbers
    // Index 2 magnetic quantum numbers
    int m2 = (is_valid_iter) ? im2 - l2po : 0;
    int ilm2 = (is_valid_iter) ? il2 + m2 : 0;
    int ilm2e = (is_valid_iter) ? ilm2 + ndi : 0;
    int i2 = (is_valid_iter) ? in2 + ilm2 : 0;
    int i2e = (is_valid_iter) ? in2 + ilm2e : 0;
    int j2 = (is_valid_iter) ? in1 + ilm2 : 0;
    int j2e = (is_valid_iter) ? in1 + ilm2e : 0;
    // End of index 2 magnetic quantum numbers
    magmaDoubleComplex cgh, cgk, zvalue;
    cgh = (is_valid_iter) ?
	magma_ghit<0>(
	  0, nbl, l1, m1, l2, m2, rxx, ryy, rzz, ind3j, v3j0, vh,
	  vyhj, vj0, vyj0, vj, li, le, rac3j
      ) :
	cz0;
    cgk = (is_valid_iter) ?
	magma_ghit<0>(
        1, nbl, l1, m1, l2, m2, rxx, ryy, rzz, ind3j, v3j0, vh,
	  vyhj, vj0, vyj0, vj, li, le, rac3j
      ) :
	cz0;
    if (is_valid_iter) {
	vec_am[(i1 - 1) * ndit + i2 - 1] = cgh;
	vec_am[(i1 - 1) * ndit + i2e - 1] = cgk;
	vec_am[(i1e - 1) * ndit + i2 - 1] = cgk;
	vec_am[(i1e - 1) * ndit + i2e - 1] = cgh;
	zvalue = magma_zprod(cgh, rsh);
	vec_am[(j1 - 1) * ndit + j2 - 1] = zvalue;
	vec_am[(j1e - 1) * ndit + j2e - 1] = zvalue;
	zvalue = magma_zprod(cgk, rsk);
	vec_am[(j1 - 1) * ndit + j2e - 1] = zvalue;
	vec_am[(j1e - 1) * ndit + j2 - 1] = zvalue;
    }
  }
  
  magma_int_t diag_iters = nsph * li * max_litpo;
#pragma omp target teams distribute parallel for \
  device(device_id)
  for(magma_int_t iter = 0; iter < diag_iters; ++iter) {
    magma_int_t t = iter;
    int im1 = (t % max_litpo) + 1;
    t /= max_litpo;
    int l1 = (t % li) + 1;
    t /= li;
    int n1 = (t % nsph) + 1;
    magmaDoubleComplex result_e, result_m;
    int in1 = (n1 - 1) * nlim;
    int rindex = (l1 - 1) * nsph + n1 - 1;
    magmaDoubleComplex dm = rmi[rindex];
    magmaDoubleComplex de = rei[rindex];
    int l1po = l1 + 1;
    int il1 = l1po * l1;
    int l1tpo = l1po + l1;
    int m1 = im1 - l1po;
    int ilm1 = il1 + m1;
    int i1 = in1 + ilm1;
    int i1e = i1 + ndi;
    result_m = dm;
    result_e = de;
    if (im1 <= l1tpo) {
	vec_am[(i1 - 1) * ndit + i1 -1] = result_m;
	vec_am[(i1e - 1) * ndit + i1e -1] = result_e;
    }
  } // iter loop
#pragma omp taskwait
  return MAGMA_SUCCESS;
}

#pragma omp begin declare target device_type(any)
template<int IHI> magmaDoubleComplex magma_ghit(
  int ipamo, int nbl, int l1, int m1, int l2, int m2,
  double *rxx, double *ryy, double *rzz, int *ind3j,
  double *v3j0, magmaDoubleComplex *vh, magmaDoubleComplex *vyhj,
  magmaDoubleComplex *vj0, magmaDoubleComplex *vyj0,
  magmaDoubleComplex vj, int li, int le, double rac3j[]
) {
  /* NBL identifies transfer vector going from N2 to N1;
   * IHI=0 for Hankel, IHI=1 for Bessel, IHI=2 for Bessel from origin;
   * depending on IHI, IPAM=0 gives H or I, IPAM= 1 gives K or L. */
  const magmaDoubleComplex cz0 = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex uim = MAGMA_Z_MAKE(0.0, 1.0);
  
  // label 10
  int lm = (le > li) ? le : li;
  int l1mp = l1 - ipamo;
  int l1po = l1 + 1;
  int m1mm2 = m1 - m2;
  int m1mm2m = (m1mm2 > 0) ? m1mm2 + 1 : 1 - m1mm2;
  int lminpo = (l2 - l1mp > 0) ? l2 - l1mp + 1 : l1mp - l2 + 1;
  int lmaxpo = l2 + l1mp + 1;
  int i3j0in = ind3j[lm * l1mp + l2 - 1];
  int ilin = (m1mm2m > lminpo && (m1mm2m - lminpo) % 2 != 0) ? 0 : -1;
  int isn = (m1 % 2 != 0) ? -1 : 1;
  if (lminpo % 2 == 0) {
    isn *= -1;
    if (l2 > l1mp) isn *= -1;
  }
  // label 12
  int nblmo = nbl - 1;
  
  if constexpr (IHI == 0) {
    int litpo = li + li + 1;
    int litpos = litpo * litpo;
    int nbhj = nblmo * litpo;
    int nby = nblmo * litpos;
    magmaDoubleComplex result = cz0;
    for (int jm24 = 1; jm24 <= 3; jm24++) {
      magmaDoubleComplex csum = cz0;
      int mu = jm24 - 2;
      int mupm1 = mu + m1;
      int mupm2 = mu + m2;
      if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	int jsn = -isn;
	if (mu == 0) jsn = isn;
	double cr = magma_cgev(ipamo, mu, l1, m1) * magma_cgev(0, mu, l2, m2);
	int i3j0 = i3j0in;
	if (mupm1 == 0 && mupm2 == 0) {
	  int lt14 = lminpo;
	  while (lt14 <= lmaxpo) {
	    i3j0++;
	    int l3 = lt14 - 1;
	    int ny = l3 * l3 + lt14;
	    double aors = 1.0 * (l3 + lt14);
	    double f3j = v3j0[i3j0 - 1] * v3j0[i3j0 - 1] * sqrt(aors) * jsn;
	    // dcomplex cfun = (vh[nbhj + lt14 - 1] * vyhj[nby + ny - 1]) * f3j;
	    magmaDoubleComplex cfun = magma_zprod(vh[nbhj + lt14 - 1], vyhj[nby + ny - 1]);
	    cfun = magma_zscale(cfun, f3j);
	    csum = magma_zadd(csum, cfun);
	    jsn *= -1;
	    lt14 += 2;
	  }
	  // goes to 22
	} else { // label 16
	  magma_r3jjr(l1mp, l2, -mupm1, mupm2, rac3j);
	  int il = ilin;
	  int lt20 = lminpo;
	  while (lt20 <= lmaxpo) {
	    i3j0++;
	    if (m1mm2m <= lt20) {
	      il += 2;
	      int l3 = lt20 - 1;
	      int ny = l3 * l3  + lt20 + m1mm2;
	      double aors = 1.0 * (l3 + lt20);
	      double f3j = rac3j[il - 1] * v3j0[i3j0 - 1] * sqrt(aors) * jsn;
	      // dcomplex cfun = (vh[nbhj + lt20 - 1] * vyhj[nby + ny - 1]) * f3j;
	      magmaDoubleComplex cfun = magma_zprod(vh[nbhj + lt20 - 1], vyhj[nby + ny - 1]);
	      cfun = magma_zscale(cfun, f3j);
	      csum = magma_zadd(csum, cfun);
	    }
	    // label 20
	    jsn *= -1;
	    lt20 += 2;
	  }
	}
	// label 22
	//csum *= cr;
	csum = magma_zscale(csum, cr);
	result = magma_zadd(result, csum);
      }
      // Otherwise there is nothing to add
    } // jm24 loop. Should go to 70
    magmaDoubleComplex cr = (ipamo == 1) ?
      MAGMA_Z_MAKE(0.0, sqrt(FOUR_PI * (l1 + l1mp) * (l1 + l1po) * (l2 + l2 + 1) / l1po)) :
      MAGMA_Z_MAKE(sqrt(FOUR_PI * (l1 + l1po) * (l2 + l2 + 1)), 0.0);
    return magma_zprod(result, cr);
  } else if constexpr (IHI == 1) {
    int litpo = li + li + 1;
    int litpos = litpo * litpo;
    int nbhj = nblmo * litpo;
    int nby = nblmo * litpos;
    magmaDoubleComplex result = cz0;
    for (int jm44 = 1; jm44 <= 3; jm44++) {
      magmaDoubleComplex csum = cz0;
      int mu = jm44 - 2;
      int mupm1 = mu + m1;
      int mupm2 = mu + m2;
      if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	int jsn = - isn;
	if (mu == 0) jsn = isn;
	double cr = magma_cgev(ipamo, mu, l1, m1) * magma_cgev(0, mu, l2, m2);
	int i3j0 = i3j0in;
	if (mupm1 == 0 && mupm2 == 0) {
	  int lt34 = lminpo;
	  while (lt34 <= lmaxpo) {
	    i3j0++;
	    int l3 = lt34 - 1;
	    int ny = l3 * l3 + lt34;
	    double aors = 1.0 * (l3 + lt34);
	    double f3j = v3j0[i3j0 - 1] * v3j0[i3j0 - 1] * sqrt(aors) * jsn;
	    magmaDoubleComplex cfun = magma_zprod(vj, vyhj[nby + ny - 1]);
	    cfun = magma_zscale(cfun, f3j);
	    csum = magma_zadd(csum, cfun);
	    jsn *= -1;
	    lt34 += 2;
	  }
	  // goes to 42
	} else { // label 36
	  magma_r3jjr(l1mp, l2, -mupm1, mupm2, rac3j);
	  int il = ilin;
	  int lt40 = lminpo;
	  while (lt40 <= lmaxpo) {
	    i3j0++;
	    if (m1mm2m <= lt40) {
	      il += 2;
	      int l3 = lt40 - 1;
	      int ny = l3 * l3  + lt40 + m1mm2;
	      double aors = 1.0 * (l3 + lt40);
	      double f3j = rac3j[il - 1] * v3j0[i3j0 - 1] * sqrt(aors) * jsn;
	      magmaDoubleComplex cfun = magma_zprod(vj, vyhj[nby + ny - 1]);
	      cfun = magma_zscale(cfun, f3j);
	      csum = magma_zadd(csum, cfun);
	    }
	    // label 40
	    jsn *= -1;
	    lt40 += 2;
	  }
	}
	// label 42
	csum = magma_zscale(csum, cr);
	result = magma_zadd(result, csum);
      }
      // Otherwise there is nothing to add
    } // jm44 loop. Should go to 70
    magmaDoubleComplex cr = (ipamo == 1) ?
      MAGMA_Z_MAKE(0.0, sqrt(FOUR_PI * (l1 + l1mp) * (l1 + l1po) * (l2 + l2 + 1) / l1po)) :
      MAGMA_Z_MAKE(sqrt(FOUR_PI * (l1 + l1po) * (l2 + l2 + 1)), 0.0);
    return magma_zprod(result, cr);
  } else { // The project grants that IHI = 2 whenever it is neither 0 nor 1.
    int lmtpo = li + le + 1;
    int lmtpos = lmtpo * lmtpo;
    int nbhj = nblmo * lmtpo;
    int nby = nblmo * lmtpos;
    magmaDoubleComplex result = cz0;
    if (rxx[nbl - 1] == 0.0 && ryy[nbl - 1] == 0.0 && rzz[nbl - 1] == 0.0) {
      if (ipamo == 0) {
	if (l1 == l2 && m1 == m2) return MAGMA_Z_MAKE(1.0, 0.0);
      }
      return MAGMA_Z_MAKE(0.0, 0.0);
    }
    for (int jm64 = 1; jm64 <= 3; jm64++) {
      magmaDoubleComplex csum = cz0;
      int mu = jm64 - 2;
      int mupm1 = mu + m1;
      int mupm2 = mu + m2;
      if (mupm1 >= -l1mp && mupm1 <= l1mp && mupm2 >= - l2 && mupm2 <= l2) {
	int jsn = -isn;
	if (mu == 0) jsn = isn;
	double cr = magma_cgev(ipamo, mu, l1, m1) * magma_cgev(0, mu, l2, m2);
	int i3j0 = i3j0in;
	if (mupm1 == 0 && mupm2 == 0) {
	  int lt54 = lminpo;
	  while (lt54 <= lmaxpo) {
	    i3j0++;
	    int l3 = lt54 - 1;
	    int ny = l3 * l3 + lt54;
	    double aors = 1.0 * (l3 + lt54);
	    double f3j = v3j0[i3j0 - 1] * v3j0[i3j0 - 1] * sqrt(aors) * jsn;
	    magmaDoubleComplex cfun = magma_zprod(vj0[nbhj + lt54 - 1], vyj0[nby + ny - 1]);
	    cfun = magma_zscale(cfun, f3j);
	    csum = magma_zadd(csum, cfun);
	    jsn *= -1;
	    lt54 += 2;
	  }
	  // goes to 62
	} else { // label 56
	  magma_r3jjr(l1mp, l2, -mupm1, mupm2, rac3j);
	  int il = ilin;
	  int lt60 = lminpo;
	  while (lt60 <= lmaxpo) {
	    i3j0++;
	    if (m1mm2m <= lt60) {
	      il += 2;
	      int l3 = lt60 - 1;
	      int ny = l3 * l3  + lt60 + m1mm2;
	      double aors = 1.0 * (l3 + lt60);
	      double f3j = rac3j[il - 1] * v3j0[i3j0 - 1] * sqrt(aors) * jsn;
	      magmaDoubleComplex cfun = magma_zprod(vj0[nbhj + lt60 - 1], vyj0[nby + ny - 1]);
	      cfun = magma_zscale(cfun, f3j);
	      csum = magma_zadd(csum, cfun);
	    }
	    // label 60
	    jsn *= -1;
	    lt60 += 2;
	  }
	}
	// label 62
	csum = magma_zscale(csum, cr);
	result = magma_zadd(result, csum);
      }
      // Otherwise there is nothing to add
    } // jm64 loop. Should go to 70
    magmaDoubleComplex cr = (ipamo == 1) ?
      MAGMA_Z_MAKE(0.0, sqrt(FOUR_PI * (l1 + l1mp) * (l1 + l1po) * (l2 + l2 + 1) / l1po)) :
      MAGMA_Z_MAKE(sqrt(FOUR_PI * (l1 + l1po) * (l2 + l2 + 1)), 0.0);
    return magma_zprod(result, cr);
  }
}
#pragma omp end declare target

#pragma omp begin declare target device_type(any)
/**
 * \brief Compute the 3j symbol for Clebsch-Gordan coefficients for JJ transitions.
 *
 * This function calculates the 3j(J,J2,J3;-M2-M3,M2,M3) symbol for the Clebsch-Gordan
 * coefficients. See Appendix a.3.1 in Borghese, Denti & Saija (2007). This function
 * is a specialized version to work on GPUs with MAGMA.
 *
 * \param j2: `int` Value of J2. [IN]
 * \param j3: `int` Value of J3. [IN]
 * \param m2: `int` Value of M2. [IN]
 * \param m3: `int` Value of M3. [IN]
 * \param rac3j: `double[]` Vector of 3j symbols. [OUT]
 */
void magma_r3jjr(int j2, int j3, int m2, int m3, double rac3j[]) {
  int jmx = j3 + j2;
  int jdf = j3 - j2;
  int m1 = -m2 - m3;
  int abs_jdf = (jdf >= 0) ? jdf : -jdf;
  int abs_m1 = (m1 >= 0) ? m1 : -m1;
  int jmn = (abs_jdf > abs_m1) ? abs_jdf : abs_m1;
  int njmo = jmx - jmn;
  int jf = jmx + jmx + 1;
  int isn = ((jdf + m1) % 2 == 0) ? 1 : -1;
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
	  rac3j[irr - 2] = -(rac3j[irr - 1] * dj + rac3j[irr] * cjp * j1) / (cj * j1po);
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
	    rac3j[irl70 - 1] * dj + rac3j[irl70 - 2] * cj * j1po
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
#pragma omp end declare target

void magma_zinvert_resident(
  magmaDoubleComplex *a, magma_int_t n, int &jer, magma_queue_t queue,
  int device_id, const RuntimeSettings& rs
) {
  string message;
  char buffer[128];
  magma_int_t err = MAGMA_SUCCESS, ref_err = MAGMA_SUCCESS;
  magma_int_t info;
  magma_int_t m = (magma_int_t)n; // changed rows; a - mxm matrix
  magma_int_t mm = m * m; // size of a
  const magmaDoubleComplex magma_zero = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex magma_one = MAGMA_Z_MAKE(1.0, 0.0);
  magmaDoubleComplex *h_a = a;
#pragma omp target data use_device_ptr(a) device(device_id)
  {
    if (rs.invert_mode == RuntimeSettings::INV_MODE_LU) { 
      // >>> LU INVERSION <<<
      magmaDoubleComplex *dwork; // work space pointer on device
      magma_int_t ldwork = m * magma_get_zgetri_nb(m); // optimal block size
      err = magma_zmalloc(&dwork, ldwork); // dev. mem. for ldwork
      magma_int_t *piv = new magma_int_t[m];
      if (err != MAGMA_SUCCESS) {
	message = "ERROR: could not allocate dwork!\n";
	rs.logger->err(message);
	exit(1);
      }
      magma_zgetrf_gpu(m, m, a, m, piv, &info);
      magma_zgetri_gpu(m, a, m, piv, dwork, ldwork, &info);
      delete[] piv; // delete piv created by magma_zgetrf()
      magma_free(dwork);
      if (rs.use_refinement) {
	// magma_newton makes a back-up copy of the unrefined inverted
	// matrix on the host. If refinement fails, the err flag is set
	// to a non-zero value, to prevent corrupting the inverted matrix.
	ref_err = magma_newton(rs, h_a, m, a, queue);
	if (ref_err != 0) {
	  err = ref_err;
	}
      }
      // >>> END OF LU INVERSION <<<
    } else if (rs.invert_mode == RuntimeSettings::INV_MODE_GESV) {
      // >>> GESV INVERSION <<<
      magmaDoubleComplex *h_id = new magmaDoubleComplex[mm]();
      magmaDoubleComplex *d_id;
      magma_int_t *piv = new magma_int_t[m];
      for (magma_int_t i = 0; i < m; i++) {
	h_id[i * (m + 1)] = magma_one;
      }
      err = magma_zmalloc(&d_id, mm);
      if (err != MAGMA_SUCCESS) {
	message = "ERROR: could not allocate d_id!\n";
	rs.logger->err(message);
	exit(1);
      }
      magma_zsetmatrix(m, m, h_id, m, d_id , m, queue); // copy identity matrix to device
      delete[] h_id;
      magma_zgesv_gpu(m, m, a, m, piv, d_id, m, &info);
      delete[] piv; // free host memory
      magma_zcopy(mm, d_id, 1, a, 1, queue); // copy d_id -> a;
      magma_free(d_id);
      if (rs.use_refinement) {
	// magma_newton makes a back-up copy of the unrefined inverted
	// matrix on the host. If refinement fails, the err flag is set
	// to a non-zero value, to prevent corrupting the inverted matrix.
	ref_err = magma_newton(rs, h_a, m, a, queue);
	if (ref_err != 0) {
	  err = ref_err;
	}
      }
      if (ref_err == MAGMA_SUCCESS) { // Refinement did improve the inversion accuracy, so we retireve the refined matrix.
	magma_zgetmatrix(m, m, a, m, h_a, m, queue); // copy d_id -> a
      }
      // >>> END OF GESV INVERSION <<<
    }
  } // end of target region
  jer = (int)err;
}

magma_int_t magma_ztm(
  magmaDoubleComplex *vec_am, ParticleDescriptor *c1, double *rxx, double *ryy,
  double *rzz, int *ind3j, double *v3j0, magmaDoubleComplex *vh, magmaDoubleComplex *vyhj,
  magmaDoubleComplex *vj0, magmaDoubleComplex *vyj0, magmaDoubleComplex *sam_v,
  magmaDoubleComplex *gis_v, magmaDoubleComplex *gls_v, magmaDoubleComplex *vec_am0m,
  int device_id
) {
  magma_int_t result = MAGMA_SUCCESS;
  const magmaDoubleComplex cz0 = MAGMA_Z_MAKE(0.0, 0.0);
  const magma_int_t k2max = c1->li * (c1->li + 2);
  const magma_int_t k3max = c1->le * (c1->le + 2);
  const magma_int_t nsph = c1->nsph;
  const magma_int_t ncou = nsph * (nsph - 1);
  const magma_int_t li = c1->li;
  const magma_int_t litpo = li + li + 1;
  const magma_int_t litpos = litpo * litpo;
  const magma_int_t nlim = li * (li + 2);
  const magma_int_t ndi = nsph * li * (li + 2);
  const magma_int_t ndit = ndi + ndi;
  const magma_int_t le = c1->le;
  const magma_int_t nlem = c1->nlem;
  const magma_int_t nlemt = nlem + nlem;
  const magma_int_t lmtpo = c1->lmtpo;
  const magma_int_t lmtpos = c1->lmtpos;
  const magma_int_t ind3j_size = (c1->lm + 1) * c1->lm;
  const magma_int_t nv3j = c1->nv3j;
  magmaDoubleComplex vj = MAGMA_Z_MAKE(real(c1->vj), imag(c1->vj));

  for (magma_int_t n2 = 1; n2 <= nsph; n2++) {
    for (magma_int_t k2 = 1; k2 <= k2max; k2++) {
      for (magma_int_t k3 = 1; k3 <= k3max; k3++) {
	double rac3j[128];
	magma_int_t l2 = (magma_int_t)sqrt(k2 + 1);
	magma_int_t im2 = k2 - (l2 * l2) + 1;
	if (im2 == 0) {
	  l2--;
	  im2 = 2 * l2+1;
	}
	else if (im2 > 2 * l2 + 1) {
	  im2 -= 2 * l2 + 1;
	  l2++;
	}
	magma_int_t l3 = (magma_int_t)sqrt(k3 + 1);
	magma_int_t im3 = k3 - (l3 * l3) + 1;
	if (im3 == 0) {
	  l3--;
	  im3 = 2 * l3 + 1;
	}
	else if (im3 > 2 * l3 + 1) {
	  im3 -= 2 * l3 + 1;
	  l3++;
	}
	magma_int_t i2 = (n2 - 1) * li * (li + 2) + l2 * l2 + im2 - 1;
	magma_int_t m2 = -l2 - 1 + im2;
	magma_int_t i3 = l3 * l3 + im3 - 1;
	magma_int_t m3 = -l3 - 1 + im3;
	magma_int_t vec_index = (i2 - 1) * nlem + i3 - 1;
	// gis_v[vec_index] = ghit(2, 0, n2, l2, m2, l3, m3, c1, rac3j);
	gis_v[vec_index] = magma_ghit<2>(
	  0, n2, l2, m2, l3, m3, rxx, ryy, rzz, ind3j,
	  v3j0, vh, vyhj, vj0, vyj0, vj, li, le, rac3j
	);
	//gls_v[vec_index] = ghit(2, 1, n2, l2, m2, l3, m3, c1, rac3j);
	gls_v[vec_index] = magma_ghit<2>(
	  1, n2, l2, m2, l3, m3, rxx, ryy, rzz, ind3j,
	  v3j0, vh, vyhj, vj0, vyj0, vj, li, le, rac3j
	);
      } // close k3 loop, former l3 + im3 loops
    } // close k2 loop, former l2 + im2 loops
  } // close n2 loop

  for (magma_int_t i1 = 1; i1 <= ndi; i1++) {
    for (magma_int_t i3 = 1; i3 <= nlem; i3++) {
      magmaDoubleComplex sum1 = cz0;
      magmaDoubleComplex sum2 = cz0;
      magmaDoubleComplex sum3 = cz0;
      magmaDoubleComplex sum4 = cz0;
      magma_int_t i1e = i1 + ndi;
      magma_int_t i3e = i3 + nlem;
      for (magma_int_t i2 = 1; i2 <= ndi; i2++) {
	magmaDoubleComplex pr;
	magma_int_t i2e = i2 + ndi;
	magma_int_t vec_ind_g23 = (i2 - 1) * nlem + i3 - 1;
	magmaDoubleComplex gie = gis_v[vec_ind_g23];
	magmaDoubleComplex gle = gls_v[vec_ind_g23];
	magma_int_t vec_ind_a1 = (i1 - 1) * ndit;
	magma_int_t vec_ind_a1e = (i1 - 1 + ndi) * ndit;
	magmaDoubleComplex a1 = vec_am[vec_ind_a1 + i2 - 1];
	magmaDoubleComplex a2 = vec_am[vec_ind_a1 + i2e - 1];
	magmaDoubleComplex a3 = vec_am[vec_ind_a1e + i2 - 1];
	magmaDoubleComplex a4 = vec_am[vec_ind_a1e + i2e - 1];
	// sum1 += (a1 * gie + a2 * gle);
	pr = magma_zprod(a1, gie);
	sum1 = magma_zadd(sum1, pr);
	pr = magma_zprod(a2, gle);
	sum1 = magma_zadd(sum1, pr);
	// sum2 += (a1 * gle + a2 * gie);
	pr = magma_zprod(a1, gle);
	sum2 = magma_zadd(sum2, pr);
	pr = magma_zprod(a2, gie);
	sum2 = magma_zadd(sum2, pr);
	// sum3 += (a3 * gie + a4 * gle);
	pr = magma_zprod(a3, gie);
	sum3 = magma_zadd(sum3, pr);
	pr = magma_zprod(a4, gle);
	sum3 = magma_zadd(sum3, pr);
	// sum4 += (a3 * gle + a4 * gie);
	pr = magma_zprod(a3, gle);
	sum4 = magma_zadd(sum4, pr);
	pr = magma_zprod(a4, gie);
	sum4 = magma_zadd(sum4, pr);
      } // i2 loop
      magma_int_t vec_ind1 = (i1 - 1) * nlemt;
      magma_int_t vec_ind1e = (i1e - 1) * nlemt;
      sam_v[vec_ind1 + i3 - 1] = sum1;
      sam_v[vec_ind1 + i3e - 1] = sum2;
      sam_v[vec_ind1e + i3 - 1] = sum3;
      sam_v[vec_ind1e + i3e - 1] = sum4;
    } // i3 loop
  } // i1 loop
  
  for (magma_int_t i1 = 1; i1 <= ndi; i1++) {
    for (magma_int_t i0 = 1; i0 <= nlem; i0++) {
      magma_int_t vec_index = (i1 - 1) * nlem + i0 - 1;
      gis_v[vec_index] = MAGMA_Z_CONJ(gis_v[vec_index]);
      gls_v[vec_index] = MAGMA_Z_CONJ(gls_v[vec_index]);
    } // i0 loop
  } // i1 loop
  
  for (magma_int_t i0 = 1; i0 <= nlem; i0++) {
    for (magma_int_t i3 = 1; i3 <= nlemt; i3++) {
      magma_int_t i0e = i0 + nlem;
      magmaDoubleComplex sum1 = cz0;
      magmaDoubleComplex sum2 = cz0;
      for (magma_int_t i1 = 1; i1 <= ndi; i1 ++) {
	magmaDoubleComplex pr;
	magma_int_t i1e = i1 + ndi;
	magma_int_t vec_ind1 = (i1 - 1) * nlemt;
	magma_int_t vec_ind1e = (i1e - 1) * nlemt;
	magmaDoubleComplex a1 = sam_v[vec_ind1 + i3 - 1];
	magmaDoubleComplex a2 = sam_v[vec_ind1e + i3 - 1];
	magma_int_t vec_index = (i1 - 1) * nlem + i0 - 1;
	magmaDoubleComplex gie = gis_v[vec_index];
	magmaDoubleComplex gle = gls_v[vec_index];
	// sum1 += (a1 * gie + a2 * gle);
	pr = magma_zprod(a1, gie);
	sum1 = magma_zadd(sum1, pr);
	pr = magma_zprod(a2, gle);
	sum1 = magma_zadd(sum1, pr);
	// sum2 += (a1 * gle + a2 * gie);
	pr = magma_zprod(a1, gle);
	sum2 = magma_zadd(sum2, pr);
	pr = magma_zprod(a2, gie);
	sum2 = magma_zadd(sum2, pr);
      } // i1 loop
      magma_int_t vec_ind0 = (i0 - 1) * nlemt;
      magma_int_t vec_ind0e = (i0e - 1) * nlemt;
      vec_am0m[vec_ind0 + i3 - 1] = MAGMA_Z_NEGATE(sum1);
      vec_am0m[vec_ind0e + i3 - 1] = MAGMA_Z_NEGATE(sum2);
    } // i3 loop
  } // i0 loop
  return result;
}

#endif // USE_TARGET_OFFLOAD

void magma_zinvert(dcomplex **mat, np_int n, int &jer, int device_id, const RuntimeSettings& rs) {
  magmaDoubleComplex *a = (magmaDoubleComplex *)mat[0]; // pointer to first element on host
  string message;
  char buffer[128];
  magma_int_t err = MAGMA_SUCCESS, ref_err = MAGMA_SUCCESS;
  magma_queue_t queue = NULL;
  magma_device_t dev = (magma_device_t)device_id;
  magma_queue_create(dev, &queue);
  magma_int_t info;
  magma_int_t m = (magma_int_t)n; // changed rows; a - mxm matrix
  magma_int_t mm = m * m; // size of a
  const magmaDoubleComplex magma_zero = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex magma_one = MAGMA_Z_MAKE(1.0, 0.0);
  if (rs.invert_mode == RuntimeSettings::INV_MODE_LU) { 
    // >>> LU INVERSION <<<
    magmaDoubleComplex *d_a; // pointer to first element on device
    magmaDoubleComplex *dwork; // work space pointer on device
    magma_int_t ldwork = m * magma_get_zgetri_nb(m); // optimal block size
    err = magma_zmalloc(&dwork, ldwork); // dev. mem. for ldwork
    magma_int_t *piv = new magma_int_t[m];
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate dwork!\n";
      rs.logger->err(message);
      exit(1);
    }
    err = magma_zmalloc(&d_a, mm); // dev. mem. for matrix
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
    magma_zgetrf_gpu(m, m, d_a, m, piv, &info);
    magma_zgetri_gpu(m, d_a, m, piv, dwork, ldwork, &info);
    delete[] piv; // delete piv created by magma_zgetrf()
    magma_free(dwork);
    if (rs.use_refinement) {
      // magma_newton makes a back-up copy of the unrefined inverted
      // matrix on the host. If refinement fails, the err flag is set
      // to a non-zero value, to prevent corrupting the inverted matrix.
      ref_err = magma_newton(rs, a, m, d_a, queue);
    }
    if (ref_err == MAGMA_SUCCESS) { // Refinement did improve the inversion accuracy, so we retireve the refined matrix.
      magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy d_a -> a
    }
    magma_free(d_a);
    // >>> END OF LU INVERSION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_GESV) {
    // >>> GESV INVERSION <<<
    magmaDoubleComplex *h_id = new magmaDoubleComplex[mm]();
    magmaDoubleComplex *d_a, *d_id;
    magma_int_t *piv = new magma_int_t[m];
    for (magma_int_t i = 0; i < m; i++) {
      h_id[i * (m + 1)] = magma_one;
    }
    err = magma_zmalloc(&d_id, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_id!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, h_id, m, d_id , m, queue); // copy identity matrix to device
    delete[] h_id;
    // allocate matrix on device
    err = magma_zmalloc(&d_a, mm); // device memory for a
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
    magma_zgesv_gpu(m, m, d_a, m, piv, d_id, m, &info);
    delete[] piv; // free host memory
    magma_free(d_a);
    if (rs.use_refinement) {
      // magma_newton makes a back-up copy of the unrefined inverted
      // matrix on the host. If refinement fails, the err flag is set
      // to a non-zero value, to prevent corrupting the inverted matrix.
      ref_err = magma_newton(rs, a, m, d_id, queue);
    }
    if (ref_err == MAGMA_SUCCESS) { // Refinement did improve the inversion accuracy, so we retireve the refined matrix.
      magma_zgetmatrix(m, m, d_id , m, a, m, queue); // copy d_id -> a
    }
    magma_free(d_id);
    // >>> END OF GESV INVERSION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_RBT) {
    // >>> RBT INVERSION <<<
    magmaDoubleComplex *d_a; // pointer to first element on device
    magma_bool_t iterate = (magma_bool_t)rs.use_refinement;
    magma_int_t *piv = new magma_int_t[m]; // host mem.
    magmaDoubleComplex *h_id = new magmaDoubleComplex[mm]();
    for (np_int i = 0; i < m; i++) {
      h_id[i * (m + 1)] = magma_one;
    }
    err = magma_zmalloc(&d_a, mm); // device memory for a
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, a, m, d_a , m, queue); // copy a -> d_a
    magma_zgesv_rbt(iterate, m, m, a, m, h_id, m, &info);
    memcpy(a, h_id, mm * sizeof(magmaDoubleComplex));
    delete[] h_id;
    delete[] piv;
    magma_free(d_a); 
    // >>> END OF RBT INVERSION <<<
  } else if (rs.invert_mode == RuntimeSettings::INV_MODE_SVD) {
    // >>> SVD INVERSION <<<
    // Step 1: compute singular value decomposition
    magma_int_t lwork;
    magmaDoubleComplex *h_work = new magmaDoubleComplex[m];
    magma_zgesdd(
      MagmaAllVec, m, m, NULL, m, NULL, NULL, m, NULL, m, h_work, -1, NULL, NULL, &info
    );
    if (info == 0) {
      lwork = (magma_int_t)(MAGMA_Z_REAL(h_work[0]));
      message = "DEBUG: lwork = " + to_string(lwork) + "\n";
      rs.logger->log(message, LOG_DEBG);
    } else {
      message = "ERROR: magma_zgesdd failed on resource querying call!\n";
      rs.logger->err(message);
      exit(1);
    }
    delete[] h_work;
    h_work = new magmaDoubleComplex[lwork];
    magmaDoubleComplex *h_a, *h_U, *h_Vt;
    h_a = new magmaDoubleComplex[mm];
    memcpy(h_a, a, mm * sizeof(magmaDoubleComplex));
    h_U = new magmaDoubleComplex[mm];
    h_Vt = new magmaDoubleComplex[mm];
    double *h_s = new double[m];
    double *h_rwork = new double[5 * mm + 7 * m];
    magma_int_t *h_iwork = new magma_int_t[8 * m];
    magma_zgesdd(
      MagmaAllVec, m, m, h_a, m, h_s, h_U, m, h_Vt, m, h_work, lwork, h_rwork, h_iwork, &info
    );
    delete[] h_work;
    delete[] h_rwork;
    delete[] h_iwork;
    delete[] h_a;
    sprintf(buffer, "DEBUG: s[0] = %.5le; s[%s] = %.5le.\n", h_s[0], to_string(m - 1).c_str(), h_s[m - 1]);
    message = buffer;
    rs.logger->log(message, LOG_DEBG);
    
    // Step 2: Upload decomposed matix to GPU
    magmaDoubleComplex *d_U, *d_Vt, *d_a;
    err = magma_zmalloc(&d_Vt, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_Vt!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, h_Vt, m, d_Vt, m, queue);
    err = magma_zmalloc(&d_U, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_U!\n";
      rs.logger->err(message);
      exit(1);
    }
    magma_zsetmatrix(m, m, h_U, m, d_U, m, queue);
    err = magma_zmalloc(&d_a, mm);
    if (err != MAGMA_SUCCESS) {
      message = "ERROR: could not allocate d_a!\n";
      rs.logger->err(message);
      exit(1);
    }

    // Step 3: compute threshold-truncated pseudo-inverse on GPU
    double eps_mach = 2.22e-18;
    double threshold = h_s[0] * m * eps_mach; 
    for (magma_int_t i = 0; i < m; i++) {
      magmaDoubleComplex inv_s = MAGMA_Z_MAKE(h_s[i] / (h_s[i] * h_s[i] + threshold * threshold), 0.0);
      magma_zscal(m, inv_s, d_Vt + i, m, queue);
    }
    delete[] h_s;
    // d_a = d_vt^H * d_u^H
    magmablas_zgemm(
      MagmaConjTrans, MagmaConjTrans, m, m, m, magma_one, d_Vt, m, d_U, m, 
      magma_zero, d_a, m, queue
    );

    // Step 4: release matrix decomposition
    magma_free(d_U);
    magma_free(d_Vt);
    delete[] h_U;
    delete[] h_Vt;

    // Step 5: refine inversion
    if (rs.use_refinement) {
      // magma_newton makes a back-up copy of the unrefined inverted
      // matrix on the host. If refinement fails, the err flag is set
      // to a non-zero value, to prevent corrupting the inverted matrix.
      ref_err = magma_newton(rs, a, m, d_a, queue);
    }

    // Step 6: get result back to host
    if (ref_err == MAGMA_SUCCESS) { // Refinement did improve the inversion accuracy, so we retireve the refined matrix.
      magma_zgetmatrix(m, m, d_a , m, a, m, queue); // copy d_a -> a
    }
    magma_free(d_a);
    // >>> END OF SVD INVERSION <<<
  }
  magma_queue_destroy(queue); // destroy queue
  jer = (int)err;
}

magma_int_t magma_newton(
  const RuntimeSettings& rs, magmaDoubleComplex* a, const magma_int_t m,
  magmaDoubleComplex* d_a, magma_queue_t queue
) {
  magma_int_t err = MAGMA_SUCCESS;
  string message;
  char buffer[128];
  const int max_ref_iters = rs.max_ref_iters;
  const magmaDoubleComplex magma_zero = MAGMA_Z_MAKE(0.0, 0.0);
  const magmaDoubleComplex magma_one = MAGMA_Z_MAKE(1.0, 0.0);
  const magmaDoubleComplex magma_mone = MAGMA_Z_MAKE(-1.0, 0.0);
  const magmaDoubleComplex magma_two = MAGMA_Z_MAKE(2.0, 0.0);
  const magma_int_t mm = m * m;
  magmaDoubleComplex *d_a_orig, *d_ax, *d_r;
  magmaDoubleComplex *h_id_diag = new magmaDoubleComplex[m];
  magmaDoubleComplex *d_id_diag;
  double oldmax = 2.0e+16, curmax = 1.0e+16;
  for (magma_int_t hi = 0; hi < m; hi++)
    h_id_diag[hi] = magma_one;
  err = magma_zmalloc(&d_id_diag, m);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_id_diag!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_zsetvector(m, h_id_diag, 1, d_id_diag, 1, queue);
  delete[] h_id_diag;
  err = magma_zmalloc(&d_a_orig, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_a_orig!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_zsetmatrix(m, m, a, m, d_a_orig, m, queue); // copy a -> d_a_orig
  magma_zgetmatrix(m, m, d_a, m, a, m, queue); // copy pre-refinement d_a -> a
  err = magma_zmalloc(&d_ax, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_ax!\n";
    rs.logger->err(message);
    exit(1);
  }
  err = magma_zmalloc(&d_r, mm);
  if (err != MAGMA_SUCCESS) {
    message = "ERROR: could not allocate d_r!\n";
    rs.logger->err(message);
    exit(1);
  }
  magma_int_t maxindex = magma_izamax(mm, d_a, 1, queue) - 1;
  magma_queue_sync(queue);
  magmaDoubleComplex magmamax = magma_mone;
  magma_zgetvector(1, d_a + maxindex, 1, &magmamax, 1, queue);
  curmax = MAGMA_Z_ABS(magmamax);
  sprintf(buffer, "INFO: largest matrix value has modulus %.5le.\n", curmax);
  message = buffer;
  rs.logger->log(message);
  for (int ri = 0; ri < max_ref_iters; ri++) {
    oldmax = curmax;
    // Compute -A*X
    magmablas_zgemm(
      MagmaNoTrans, MagmaNoTrans, m, m, m, magma_mone, d_a, m,
      d_a_orig, m, magma_zero, d_ax, m, queue
    );
    // Transform -A*X into (I - A*X)
    magma_zaxpy(m, magma_one, d_id_diag, 1, d_ax, m + 1, queue);
    maxindex = magma_izamax(mm, d_ax, 1, queue) - 1;
    magma_queue_sync(queue);
    magma_zgetvector(1, d_ax + maxindex, 1, &magmamax, 1, queue);
    curmax = MAGMA_Z_ABS(magmamax);
    sprintf(buffer, "DEBUG: iteration %d has residue %.5le; target residue is %.5le.\n", ri, curmax, rs.accuracy_goal);
    message = buffer;
    rs.logger->log(message, LOG_DEBG);
    if (curmax < 0.5) { // Safe conditions for Newton-Schultz iteration.
      if (curmax < 0.95 * oldmax) { // Newton-Schultz iteration is improving and can proceed.
	// Compute R = X*(A*X - I)
	magmablas_zgemm(
          MagmaNoTrans, MagmaNoTrans, m, m, m, magma_one, d_ax, m,
	  d_a, m, magma_zero, d_r, m, queue
        );
	// Set X = X + R
	magma_zaxpy(mm, magma_one, d_r, 1, d_a, 1, queue);
	if (curmax < rs.accuracy_goal) {
	  message = "DEBUG: good news - optimal convergence achieved. Stopping.\n";
	  rs.logger->log(message, LOG_DEBG);
	  err = MAGMA_SUCCESS;
	  break; // ri for
        }
      } else {
	if (curmax > 0.1) {
	  sprintf(buffer, "INFO: iteration %d achieved limiting residue %.5le. Cannot reach goal. Reverting.\n", ri, curmax);
	  message = buffer;
	  rs.logger->log(message);
	  err = 1;
	} else {
	  sprintf(buffer, "WARN: iteration %d achieved limiting residue %.5le. Stopping.\n", ri, curmax);
	  message = buffer;
	  rs.logger->log(message, LOG_WARN);
	}
	break; // ri for
      }
    } else { // curmax > 0.5. Newton-Schultz iteration is dangerous.
      if (curmax < oldmax) {
	sprintf(buffer, "WARN: iteration %d has residue %.5le. Iterating is dangerous. Stopping.\n", ri, curmax);
	message = buffer;
	rs.logger->log(message, LOG_WARN);
      } else {
	err = 1;
	sprintf(buffer, "INFO: iteration %d has residue %.5le. Reverting to unrefined and stopping.\n", ri, curmax);
	message = buffer;
	rs.logger->log(message);
      }
      break; // ri for
    }
  } // end of ri for
  magma_free(d_a_orig);
  magma_free(d_ax);
  magma_free(d_r);
  magma_free(d_id_diag);
  return err;
}

#endif
