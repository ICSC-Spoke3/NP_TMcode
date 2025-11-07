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

/*! \file clu_subs.h
 *
 * \brief C++ porting of CLU functions and subroutines.
 *
 * This library includes a collection of functions that are used to solve the
 * scattering problem in the case of a cluster of spheres. The functions that
 * were generalized from the case of the single sphere are imported the `sph_subs.h`
 * library. As it occurs with the single sphere case functions, in most cases, the
 * results of calculations do not fall back to fundamental data types. They are
 * rather multi-component structures. In order to manage access to such variety
 * of return values, most functions are declared as `void` and they operate on
 * output arguments passed by reference.
 */

#ifndef INCLUDE_CLU_SUBS_H_
#define INCLUDE_CLU_SUBS_H_

/*! \brief Compute the asymmetry-corrected scattering cross-section of a cluster.
 *
 * This function computes the product between the geometrical asymmetry parameter and
 * the scattering cross-section, like `aps()`, but for a cluster of spheres. See Eq. (3.16)
 * in Borghese, Denti & Saija (2007).
 *
 * \param zpv: `double ****`
 * \param le: `int`
 * \param am0m: `complex double **`
 * \param w: `complex double **`
 * \param sqk: `double`
 * \param gapr: `double **`
 * \param gapp: `complex double **`
 */
void apc(
	 double ****zpv, int le, dcomplex **am0m, dcomplex **w,
	 double sqk, double **gapr, dcomplex **gapp
);

/*! \brief Compute the asymmetry-corrected scattering cross-section under random average
 * conditions.
 *
 * This function computes the product between the geometrical asymmetry parameter and
 * the scattering cross-section of a cluster using the random average directions.
 *
 * \param zpv: `double ****`
 * \param le: `int`
 * \param am0m: `complex double **`
 * \param inpol: `int` Polarization type.
 * \param sqk: `double`
 * \param gaprm: `double **`
 * \param gappm: `complex double **`
 */
void apcra(
	   double ****zpv, const int le, dcomplex **am0m, int inpol, double sqk,
	   double **gaprm, dcomplex **gappm
);

/*! \brief Complex inner product.
 *
 * This function performs the complex inner product. It is used by `lucin()`.
 *
 * \param z: `complex double`
 * \param am: `complex double **`
 * \param i: `int`
 * \param jf: `int`
 * \param k: `int`
 * \param nj: `int`
 * \param istep: `np_int` Size of rows in the matrix.
 */
dcomplex cdtp(dcomplex z, dcomplex *vec_am, int i, int jf, int k, int nj, np_int istep);

/*! \brief C++ porting of CGEV. QUESTION: description?
 *
 * \param ipamo: `int`
 * \param mu: `int`
 * \param l: `int`
 * \param m: `int`
 * \return result: `double`
 */
double cgev(int ipamo, int mu, int l, int m);

/*! \brief Build the multi-centered M-matrix of the cluster.
 *
 * This function constructs the multi-centered M-matrix of the cluster, according
 * to Eq. (5.28) of Borghese, Denti & Saija (2007).
 *
 * \param am: `complex double **`
 * \param c1: `ParticleDescriptor *`
 */
void cms(dcomplex **am, ParticleDescriptor *c1);

/*! \brief Compute orientation-averaged scattered field intensity.
 *
 * This function computes the intensity of the scattered field for the cluster,
 * averaged on the orientations. It is invoked for IAVM=1 (geometry referred to
 * the meridional plane). QUESTION: correct?
 *
 * \param vk: `double` Wave number.
 * \param exri: `double` External medium refractive index.
 * \param c1: `ParticleDescriptor *`
 */
void crsm1(double vk, double exri, ParticleDescriptor *c1);

/*! \brief Compute the transfer vector from N2 to N1.
 *
 * This function computes the transfer vector going from N2 to N1, using either
 * Hankel, Bessel or Bessel from origin functions.
 *
 * \param ihi: `int`
 * \param ipamo: `int`
 * \param nbl: `int`
 * \param l1: `int`
 * \param m1: `int`
 * \param l2: `int`
 * \param m2: `int`
 * \param c1: `ParticleDescriptor *`
 * \param rac3j: `dcomplex *`
 */
dcomplex ghit_d(
	      int ihi, int ipamo, int nbl, int l1, int m1, int l2, int m2,
	      ParticleDescriptor *c1, double *rac3j
);

/*! \brief Compute the transfer vector from N2 to N1.
 *
 * This function computes the transfer vector going from N2 to N1, using either
 * Hankel, Bessel or Bessel from origin functions.
 *
 * \param ihi: `int`
 * \param ipamo: `int`
 * \param nbl: `int`
 * \param l1: `int`
 * \param m1: `int`
 * \param l2: `int`
 * \param m2: `int`
 * \param c1: `ParticleDescriptor *` Poiunter to a ParticleDescriptor instance.
 */
dcomplex ghit(
	      int ihi, int ipamo, int nbl, int l1, int m1, int l2, int m2,
	      ParticleDescriptor *c1
);

/*! \brief Compute Hankel funtion and Bessel functions.
 *
 * This function constructs the Hankel function and the Bessel functions vectors. See
 * page 331 in Borghese, Denti & Saija (2007).
 *
 * \param exri: `double` External medium refractive index.
 * \param vk: `double` Wave number.
 * \param jer: `int &` Reference to error code flag.
 * \param lcalc: `int &` Reference to the highest order accounted for in calculation.
 * \param arg: `complex\<double\> &`
 * \param c1: `ParticleDescriptor *`
 */
void hjv(
	 double exri, double vk, int &jer, int &lcalc, dcomplex &arg,
	 ParticleDescriptor *c1
);

/*! \brief Invert the multi-centered M-matrix.
 *
 * This function performs the inversion of the multi-centered M-matrix through
 * LU decomposition. See Eq. (5.29) in Borghese, Denti & Saija (2007).
 *
 * \param am: `complex double **`
 * \param nddmst: `const int64_t`
 * \param n: `int64_t`
 * \param ier: `int &`
 */
void lucin(dcomplex **am, const np_int nddmst, np_int n, int &ier);

/*! \brief Compute the average extinction cross-section.
 *
 * This funbction computes the average extinction cross-section starting
 * from the definition of the scattering amplitude. See Sec. 3.2.1 of Borghese,
 * Denti & Saija (2007).
 *
 * \param vk: `double` Wave number.
 * \param exri: `double` External medium refractive index.
 * \param fsac: Matrix of complex
 * \param cextlr: `double **`
 * \param cext: `double **`
 */
void mextc(double vk, double exri, dcomplex **fsac, double **cextlr, double **cext);

/*! \brief Compute cross-sections and forward scattering amplitude for the cluster.
 *
 * This function computes the scattering, absorption and extinction cross-sections
 * of the cluster, together with the Forward Scattering Amplitude.
 *
 * This function is intended to evaluate the particle cross-section. QUESTIUON: correct?
 * \param vk: `double` Wave number.
 * \param exri: `double` External medium refractive index.
 * \param c1: `ParticleDescriptor *`
 */
void pcros(double vk, double exri, ParticleDescriptor *c1);

/*! \brief Compute orientation-averaged cross-sections and forward scattering amplitude.
 *
 * This function computes the orientation-averaged scattering, absorption and extinction
 * cross-sections of the cluster, together with the averaged Forward Scattering Amplitude.
 *
 * \param vk: `double` Wave number.
 * \param exri: `double` External medium refractive index.
 * \param inpol: `int` Incident field polarization type.
 * \param c1: `ParticleDescriptor *`
 */
void pcrsm0(double vk, double exri, int inpol, ParticleDescriptor *c1);

/*! \brief Transform Cartesian quantities to spherical coordinates ones.
 *
 * This function performs a conversion from the Cartesian coordinates system to the
 * spherical one. It is used by `sphar()`.
 *
 * \param x: `double` X-axis Cartesian coordinate.
 * \param y: `double` Y-axis Cartesian coordinate.
 * \param z: `double` Z-axis Cartesian coordinate.
 * \param r: `double &` Reference to radial vector (output value).
 * \param cth: `double &` Reference to the cosine of the azimuth coordinate (output value).
 * \param sth: `double &` Reference to the sine of the azimuth coordinate (output value).
 * \param cph: `double &` Reference to the cosine of the elevation coordinate (output value).
 * \param sph: `double &` Reference to the sine of the elevation coordinate (output value).
 */
void polar(
	   double x, double y, double z, double &r, double &cth, double &sth,
	   double &cph, double &sph
);

/*! \brief Compute the 3j symbol for Clebsch-Gordan coefficients towards J=0.
 *
 * This function calculates the 3j(J,J2,J3;0,0,0) symbol for the Clebsch-Gordan
 * coefficients. See Appendix a.3.1 in Borghese, Denti & Saija (2007).
 *
 * \param j2: `int`
 * \param j3: `int`
 * \param rac3j: `double *` Vector of 3j symbols.
 */
void r3j000(int j2, int j3, double *rac3j);

/*! \brief Compute the 3j symbol for Clebsch-Gordan coefficients for JJ transitions.
 *
 * This function calculates the 3j(J,J2,J3;-M2-M3,M2,M3) symbol for the Clebsch-Gordan
 * coefficients. See Appendix a.3.1 in Borghese, Denti & Saija (2007).
 *
 * \param j2: `int`
 * \param j3: `int`
 * \param m2: `int`
 * \param m3: `int`
 * \param rac3j: `double *` Vector of 3j symbols.
 */
void r3jjr(int j2, int j3, int m2, int m3, double *rac3j);

/*! \brief Compute the 3j symbol for Clebsch-Gordan coefficients for JJ transitions.
 *
 * This function calculates the 3j(J,J2,J3;-M2-M3,M2,M3) symbol for the Clebsch-Gordan
 * coefficients. See Appendix a.3.1 in Borghese, Denti & Saija (2007).
 *
 * \param j2: `int`
 * \param j3: `int`
 * \param m2: `int`
 * \param m3: `int`
 * \param rac3j: `double *`
 */
void r3jjr_d(int j2, int j3, int m2, int m3, double *rac3j);

/*! \brief Compute the 3j symbol for Clebsch-Gordan coefficients for JM transitions.
 *
 * This function calculates the 3j(J,J2,J3;M1,M,-M1-M) symbol for the Clebsch-Gordan
 * coefficients. See Appendix a.3.1 in Borghese, Denti & Saija (2007).
 *
 * \param j1: `int`
 * \param j2: `int`
 * \param j3: `int`
 * \param m1: `int`
 * \param rac3j: `double *` Vector of 3j symbols.
 */
void r3jmr(int j1, int j2, int j3, int m1, double *rac3j);

/*! \brief Compute radiation torques on a particle in Cartesian coordinates.
 *
 * This function computes radiation torque on on a cluster of spheres as the
 * result of the difference between the extinction and the scattering
 * contributions for a Cartesian coordinate system, as `rabas()`. See Sec. 4.9
 * in Borghese, Denti & Saija (2007).
 *
 * \param le: `int`
 * \param am0m: `complex double **`
 * \param w: `complex double **`
 * \param tqce: `double **`
 * \param tqcpe: `complex double **`
 * \param tqcs: `double **`
 * \param tqcps: `complex double **`
 */
void raba(
	  int le, dcomplex **am0m, dcomplex **w, double **tqce,
	  dcomplex **tqcpe, double **tqcs, dcomplex **tqcps
);

/*! \brief Compute the radiation force Cartesian components.
 *
 * This function computes the Cartesian components of the radiation force
 * exerted on a particle. See Sec. 3.2.1 in Borghese, Denti & Saija (2007).
 *
 * \param u: `double *`
 * \param up: `double *`
 * \param un: `double *`
 * \param gapv: `double *`
 * \param extins: `double`
 * \param scatts: `double`
 * \param rapr: `double &`
 * \param cosav: `double &`
 * \param fp: `double &`
 * \param fn: `double &`
 * \param fk: `double &`
 * \param fx: `double &`
 * \param fy: `double &`
 * \param fz: `double &`
 */
void rftr(
	  double *u, double *up, double *un, double *gapv, double extins, double scatts,
	  double &rapr, double &cosav, double &fp, double &fn, double &fk, double &fx,
	  double &fy, double &fz
);

/*! \brief Compute Mie cross-sections for the sphere units in the cluster.
 *
 * This function computes the scattering, absorption and extinction cross-sections
 * for the spheres composing the cluster, in terms of Mie coefficients, together
 * with the Forward Scattering Amplitude. See Sec. 4.2.1 in Borghese, Denti & Saija
 * (2007).
 *
 * \param vk: `double` Wave number
 * \param exri: `double` External medium refractive index.
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void scr0(double vk, double exri, ParticleDescriptor *c1);

/*! \brief Compute the scattering amplitude for a single sphere in an aggregate.
 *
 * This function computes the scattering amplitude for single spheres constituting
 * an aggregate. See Sec. 4.2.1 in Borghese, Denti & Saija (2007).
 *
 * \param vk: `double` Wave number.
 * \param vkarg: `double` QUESTION: definition?
 * \param exri: `double` External medium refractive index.
 * \param duk: `double *` QUESTION: definition?
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void scr2(
	  double vk, double vkarg, double exri, double *duk, ParticleDescriptor *c1
);

/*! \brief Transform sphere Cartesian coordinates to spherical coordinates.
 *
 * This function transforms the Cartesian coordinates of the spheres in an aggregate
 * to radial coordinates, then it calls `sphar()` to calculate the vector of spherical
 * harmonics of the incident field.
 *
 * \param rcf: `double **` Matrix of sphere configuration fractional radii.
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void str(double **rcf, ParticleDescriptor *c1);

/*! \brief Compute radiation torques on particles on a k-vector oriented system.
 *
 * This function computes the radiation torques resulting from the difference
 * between absorption and scattering contributions, like `rabas()`, but for a
 * coordinate system oriented along the wave vector and its orthogonal axis. See
 * Sec. 4.9 in Borghese, Denti & Saija (2007).
 *
 * \param u: `double *`
 * \param up: `double *`
 * \param un: `double *`
 * \param tqev: `double *`
 * \param tqsv: `double *`
 * \param tep: `double &`
 * \param ten: `double &`
 * \param tek: `double &`
 * \param tsp: `double &`
 * \param tsn: `double &`
 * \param tsk: `double &`
 */
void tqr(
	 double *u, double *up, double *un, double *tqev, double *tqsv, double &tep,
	 double &ten, double &tek, double &tsp, double &tsn, double &tsk
);

/*! \brief Calculate the single-centered inversion of the M-matrix.
 *
 * This function computes the single-centered inverted M-matrix appearing in Eq. (5.28)
 * of Borghese, Denti & Saija (2007).
 *
 * \param am: `complex double **`
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void ztm(dcomplex **am, ParticleDescriptor *c1);

#endif
