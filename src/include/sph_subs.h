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

/*! \file sph_subs.h
 *
 * \brief C++ porting of SPH functions and subroutines.
 *
 * This library includes a collection of functions that are used to solve the
 * scattering problem in the case of a single sphere. Some of these functions
 * are also generalized to the case of clusters of spheres. In most cases, the
 * results of calculations do not fall back to fundamental data types. They are
 * rather multi-component structures. In order to manage access to such variety
 * of return values, most functions are declared as `void` and they operate on
 * output arguments passed by reference.
 */

#ifndef INCLUDE_SPH_SUBS_H_
#define INCLUDE_SPH_SUBS_H_

/*! \brief Compute the asymmetry-corrected scattering cross-section of a sphere.
 *
 * This function computes the product between the geometrical asymmetry parameter and
 * the scattering cross-section. See Sec. 3.2.1 of Borghese, Denti & Saija (2007).
 *
 * \param zpv: `double ****` Geometrical asymmetry parameter coefficients matrix.
 * \param li: `int` Maximum field expansion order.
 * \param nsph: `int` Number of spheres.
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 * \param sqk: `double`
 * \param gaps: `double *` Geometrical asymmetry parameter-corrected cross-section.
 */
void aps(double ****zpv, int li, int nsph, ParticleDescriptor *c1, double sqk, double *gaps);

/*! \brief Complex Bessel Function.
 *
 * This function computes the complex spherical Bessel funtions \f$j\f$. It uses the
 * auxiliary functions `msta1()` and `msta2()` to determine the starting point for
 * backward recurrence. This is the `CSPHJ` implementation of the `specfun` library.
 *
 * \param n: `int` Order of the function.
 * \param z: `complex double` Argument of the function.
 * \param nm: `int &` Highest computed order.
 * \param csj: `complex double *` The desired function \f$j\f$.
 */
void cbf(int n, dcomplex z, int &nm, dcomplex *csj);

/*! \brief Clebsch-Gordan coefficients.
 *
 * This function comutes the Clebsch-Gordan coefficients for the irreducible spherical
 * tensors. See Sec. 1.6.3 and Table 1.1 in Borghese, Denti & Saija (2007).
 *
 * \param lmpml: `int`
 * \param mu: `int`
 * \param l: `int`
 * \param m: `int`
 * \return result: `double` Clebsh-Gordan coefficient.
 */
double cg1(int lmpml, int mu, int l, int m);

/*! \brief Compute the continuous variation of the refractive index and of its radial derivative.
 *
 * This function implements the continuous variation of the refractive index and of its radial
 * derivative through the materials that constitute the sphere and its surrounding medium. See
 * Sec. 5.2 in Borghese, Denti & Saija (2007).
 *
 * \param npntmo: `int`
 * \param ns: `int`
 * \param i: `int`
 * \param ic: `int`
 * \param vk: `double`
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 */
void diel(int npntmo, int ns, int i, int ic, double vk, ParticleDescriptor *c1);

/*! \brief Compute Mie scattering coefficients.
 *
 * This function determines the L-dependent Mie scattering coefficients \f$a_l\f$ and \f$b_l\f$
 * for the cases of homogeneous spheres, radially non-homogeneous spheres and, in case of sphere
 * with dielectric function, sustaining longitudinal waves. See Sec. 5.1 in Borghese, Denti
 * & Saija (2007).
 *
 * \param li: `int` Maximum field expansion order.
 * \param i: `int`
 * \param npnt: `int`
 * \param npntts: `int`
 * \param vk: `double` Wave number in scale units.
 * \param exdc: `double` External medium dielectric constant.
 * \param exri: `double` External medium refractive index.
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 * \param jer: `int &` Reference to integer error code variable.
 * \param lcalc: `int &` Reference to integer variable recording the maximum expansion order accounted for.
 * \param arg: `complex double &`
 * \param last_conf: `int` Last sphere configuration (used by CLUSTER)
 */
void dme(
	 int li, int i, int npnt, int npntts, double vk, double exdc, double exri,
	 ParticleDescriptor *c1, int &jer, int &lcalc, dcomplex &arg, int last_conf=0
);

/*! \brief Bessel function calculation control parameters.
 *
 * This function determines the control parameters required to calculate the Bessel
 * functions up to the desired precision.
 *
 * \param n: `int`
 * \param x: `double`
 * \return result: `double`
 */
double envj(int n, double x);

/*! \brief Compute the Mueller Transformation Matrix.
 *
 * This function computes the Mueller Transformation Matrix, or Phase Matrix. See
 * Sec. 2.8.1 of Borghese, Denti & Saija (2007).
 *
 * \param vint: `complex double *`
 * \param cmullr: `double **`
 * \param cmul: `double **`
 */
void mmulc(dcomplex *vint, double **cmullr, double **cmul);

/*! \brief Starting point for Bessel function magnitude.
 *
 * This function determines the starting point for backward recurrence such that
 * the magnitude of all \f$J\f$ and \f$j\f$ functions is of the order of \f$10^{-mp}\f$.
 *
 * \param x: `double` Absolute value of the argumetn to \f$J\f$ or \f$j\f$.
 * \param mp: `int` Requested order of magnitude.
 * \return result: `int` The necessary starting point.
 */
int msta1(double x, int mp);

/*! \brief Starting point for Bessel function precision.
 *
 * This function determines the starting point for backward recurrence such that
 * all \f$J\f$ and \f$j\f$ functions have `mp` significant digits.
 *
 * \param x: `double` Absolute value of the argumetn to \f$J\f$ or \f$j\f$.
 * \param n: `int` Order of the function.
 * \param mp: `int` Requested number of significant digits.
 * \return result: `int` The necessary starting point.
 */
int msta2(double x, int n, int mp);

/*! \brief Compute the amplitude of the orthogonal unit vector.
 *
 * This function computes the amplitude of the orthogonal unit vector for a geometry
 * based either on the scattering plane or on the meridional plane. It is used by `upvsp()`
 * and `upvmp()`. See Sec. 2.7 in Borghese, Denti & Saija (2007).
 *
 * \param u1: `double *`
 * \param u2: `double *`
 * \param u3: `double *`
 * \param iorth: `int`
 * \param torth: `double`
 */
void orunve(double *u1, double *u2, double *u3, int iorth, double torth);

/*! \brief Compute incident and scattered field amplitudes.
 *
 * This function computes the amplitudes of the incident and scattered field on the
 * basis of the multi-polar expansion. See Sec. 4.1.1 in Borghese, Denti and Saija (2007).
 *
 * \param up: `double *`
 * \param un: `double *`
 * \param ylm: `complex double *` Field polar spherical harmonics.
 * \param inpol: `int` Incident field polarization type (0 - linear; 1 - circular).
 * \param lw: `int`
 * \param isq: `int`
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 */
void pwma(
	  double *up, double *un, dcomplex *ylm, int inpol, int lw,
	  int isq, ParticleDescriptor *c1
);

/*! \brief Compute radiation torques on a single sphere in Cartesian coordinates.
 *
 * This function computes radiation torque on a sphere unit as the result
 * of the difference between the extinction and the scattering contributions.
 * See Sec. 4.9 in Borghese, Denti & Saija (2007).
 *
 * \param inpol: `int` Incident polarization type (0 - linear; 1 - circular)
 * \param li: `int` Maximum field expansion order.
 * \param nsph: `int` Number of spheres.
 * \param c1: `ParticleDescriptor *` Pointer to `ParticleDescriptor` data structure.
 * \param tqse: `double **`
 * \param tqspe: `complex double **`
 * \param tqss: `double **`
 * \param tqsps: `complex double **`
 */
void rabas(
	   int inpol, int li, int nsph, ParticleDescriptor *c1, double **tqse, dcomplex **tqspe,
	   double **tqss, dcomplex **tqsps
);

/*! \brief Real Bessel Function.
 *
 * This function computes the real spherical Bessel funtions \f$j\f$. It uses the
 * auxiliary functions `msta1()` and `msta2()` to determine the starting point for
 * backward recurrence. This is the `SPHJ` implementation of the `specfun` library.
 *
 * \param n: `int` Order of the function.
 * \param x: `double` Argument of the function.
 * \param nm: `int &` Highest computed order.
 * \param sj: `double[]` The desired function \f$j\f$.
 */
void rbf(int n, double x, int &nm, double sj[]);

/*! \brief Soft layer radial function and derivative.
 *
 * This function determines the radial function and its derivative for a soft layer
 * in a radially non homogeneous sphere. See Sec. 5.1 in Borghese, Denti & Saija (2007).
 *
 * \param npntmo: `int`
 * \param step: `double`
 * \param dcc: `complex double`
 * \param x: `double &`
 * \param lpo: `int`
 * \param y1: `complex double &`
 * \param y2: `complex double &`
 * \param dy1: `complex double &`
 * \param dy2: `complex double &`
 */
void rkc(
	 int npntmo, double step, dcomplex dcc, double &x, int lpo,
	 dcomplex &y1, dcomplex &y2, dcomplex &dy1, dcomplex &dy2
);

/*! \brief Transition layer radial function and derivative.
 *
 * This function determines the radial function and its derivative for a transition layer
 * in a radially non homogeneous sphere. See Sec. 5.1 in Borghese, Denti & Saija (2007).
 *
 * \param npntmo: `int`
 * \param step: `double`
 * \param x: `double &`
 * \param lpo: `int`
 * \param y1: `complex double &`
 * \param y2: `complex double &`
 * \param dy1: `complex double &`
 * \param dy2: `complex double &`
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void rkt(
	 int npntmo, double step, double &x, int lpo, dcomplex &y1,
	 dcomplex &y2, dcomplex &dy1, dcomplex &dy2, ParticleDescriptor *c1
);

/*! \brief Spherical Bessel functions.
 *
 * This function computes the spherical Bessel functions \f$y\f$. It adopts the `SPHJY`
 * implementation of the `specfun` library.
 *
 * \param n: `int` Order of the function (from 0 up).
 * \param x: `double` Argumento of the function (\f$x > 0\f$).
 * \param nm: `int &` Highest computed order.
 * \param sy: `double *` The desired function \f$y\f$.
 */
void rnf(int n, double x, int &nm, double *sy);

/*! \brief Spherical harmonics for given direction.
 *
 * This function computes the field spherical harmonics for a given direction. See Sec.
 * 1.5.2 in Borghese, Denti & Saija (2007).
 *
 * \param cosrth: `double` Cosine of direction's elevation.
 * \param sinrth: `double` Sine of direction's elevation.
 * \param cosrph: `double` Cosine of direction's azimuth.
 * \param sinrph: `double` Sine of direction's azimuth.
 * \param ll: `int` L value expansion order.
 * \param ylm: `complex double *` The requested spherical harmonics.
 */
void sphar(
	   double cosrth, double sinrth, double cosrph, double sinrph,
	   int ll, dcomplex *ylm
);

/*! \brief Compute scattering, absorption and extinction cross-sections.
 *
 * This function computes the scattering, absorption and extinction cross-sections in terms
 * of Forward Scattering Amplitudes. See Sec. 4.2.1 in Borghese, Denti & Saija (2007).
 *
 * \param tfsas: `complex double &`
 * \param nsph: `int` Number of spheres.
 * \param lm: `int` Maximum field expansion order.
 * \param vk: `double` Wave number in scale units.
 * \param exri: `double` External medium refractive index.
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 */
void sscr0(dcomplex &tfsas, int nsph, int lm, double vk, double exri, ParticleDescriptor *c1);

/*! \brief C++ Compute the scattering amplitude and the scattered field intensity.
 *
 * The role of this function is to compute the scattering amplitude and the intensity of
 * the scattered field. See Sec. 4.2 in Borghese, Denti & Saija (2007).
 *
 * \param nsph: `int` Number of spheres.
 * \param lm: `int` Maximum field expansion order.
 * \param vk: `double` Wave number in scale units.
 * \param exri: `double` External medium refractive index.
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 */
void sscr2(int nsph, int lm, double vk, double exri, ParticleDescriptor *c1);

/*! \brief Determine the geometrical asymmetry parameter coefficients.
 *
 * This function computes the coefficients that enter the definition of the geometrical
 * asymmetry parameter based on the L-value of the field expansion order. See Sec. 3.2.1
 * in Borghese, Denti & Saija (2007).
 *
 * \param lm: `int` Maximum field expansion order.
 * \param zpv: `double ****` Matrix of geometrical asymmetry parameter coefficients.
 */
void thdps(int lm, double ****zpv);

/*! \brief Compute the unitary vectors onb the polarization plane and its orthogonal
 * direction.
 *
 * This function computes the unitary vectors lying on the polarization plane and on
 * its orthogonal direction, to optimize the identification of the scattering geometry.
 * See Sec. 2.3 in Borghese, Denti & Saija (2007).
 *
 * \param thd: `double`
 * \param phd: `double`
 * \param icspnv: `int`
 * \param cost: `double`
 * \param sint: `double`
 * \param cosp: `double`
 * \param sinp: `double`
 * \param u: `double *`
 * \param up: `double *`
 * \param un: `double *`
 */
void upvmp(
	   double thd, double phd, int icspnv, double &cost, double &sint,
	   double &cosp, double &sinp, double *u, double *up, double *un
);

/*! \brief Compute the unitary vector perpendicular to incident and scattering plane.
 *
 * This function computes the unitary vector perpendicular to the incident and scattering
 * plane in a geometry based on the scattering plane. It uses `orunve()`. See Sec. 2.7 in
 * Borghese, Denti & Saija (2007).
 *
 * \param u: `double *`
 * \param upmp: `double *`
 * \param unmp: `double *`
 * \param us: `double *`
 * \param upsmp: `double *`
 * \param unsmp: `double *`
 * \param up: `double *`
 * \param un: `double *`
 * \param ups: `double *`
 * \param uns: `double *`
 * \param duk: `double *`
 * \param isq: `int &`
 * \param ibf: `int &`
 * \param scand: `double &` Scattering angle in degrees.
 * \param cfmp: `double &`
 * \param sfmp: `double &`
 * \param cfsp: `double &`
 * \param sfsp: `double &`
 */
void upvsp(
	   double *u, double *upmp, double *unmp, double *us, double *upsmp, double *unsmp,
	   double *up, double *un, double *ups, double *uns, double *duk, int &isq,
	   int &ibf, double &scand, double &cfmp, double &sfmp, double &cfsp, double &sfsp
);

/*! \brief Compute meridional plane-referred geometrical asymmetry parameter coefficients.
 *
 * This function computes the coeffcients that define the geometrical symmetry parameter
 * as defined with respect to the meridional plane. It makes use of `sphar()` and `pwma()`.
 * See Sec. 3.2.1 in Borghese, Denti & Saija (2007).
 *
 * \param iis: `int`
 * \param cost: `double` Cosine of the elevation angle.
 * \param sint: `double` Sine of the elevation angle.
 * \param cosp: `double` Cosine of the azimuth angle.
 * \param sinp: `double` Sine of the azimuth angle.
 * \param inpol: `int` Incident field polarization type (0 - linear; 1 - circular).
 * \param lm: `int` Maximum field expansion orde.
 * \param idot: `int`
 * \param nsph: `int` Number of spheres.
 * \param arg: `double *`
 * \param u: `double *`
 * \param up: `double *`
 * \param un: `double *`
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 */
void wmamp(
	   int iis, double cost, double sint, double cosp, double sinp, int inpol,
	   int lm, int idot, int nsph, double *arg, double *u, double *up,
	   double *un, ParticleDescriptor *c1
);

/*! \brief Compute the scattering plane-referred geometrical asymmetry parameter coefficients.
 *
 * This function computes the coefficients that define the geometrical asymmetry parameter based
 * on the L-value with respect to the scattering plane. It uses `sphar()` and `pwma()`. See Sec.
 * 3.2.1 in Borghese, Denti and Saija (2007).
 *
 * \param cost: `double` Cosine of elevation angle.
 * \param sint: `double` Sine of elevation angle.
 * \param cosp: `double` Cosine of azimuth angle.
 * \param sinp: `double` Sine of azimuth angle.
 * \param costs: `double` Cosine of scattering elevation angle.
 * \param sints: `double` Sine of scattering elevation angle.
 * \param cosps: `double` Cosine of scattering azimuth angle.
 * \param sinps: `double` Sine of scattering azimuth angle.
 * \param u: `double *`
 * \param up: `double *`
 * \param un: `double *`
 * \param us: `double *`
 * \param ups: `double *`
 * \param uns: `double *`
 * \param isq: `int`
 * \param ibf: `int`
 * \param inpol: `int` Incident field polarization (0 - linear; 1 -circular).
 * \param lm: `int` Maximum field expansion order.
 * \param idot: `int`
 * \param nsph: `int` Number opf spheres.
 * \param argi: `double *`
 * \param args: `double *`
 * \param c1: `ParticleDescriptor *` Pointer to a `ParticleDescriptor` data structure.
 */
void wmasp(
	   double cost, double sint, double cosp, double sinp, double costs, double sints,
	   double cosps, double sinps, double *u, double *up, double *un, double *us,
	   double *ups, double *uns, int isq, int ibf, int inpol, int lm, int idot,
	   int nsph, double *argi, double *args, ParticleDescriptor *c1
);

#endif /* SRC_INCLUDE_SPH_SUBS_H_ */
