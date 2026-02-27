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

/*! \file tra_subs.h
 *
 * \brief C++ porting of TRAPPING functions and subroutines.
 *
 * This library includes a collection of functions that are used to solve the
 * trapping problem. The functions that were generalized from the case of the
 * single sphere are imported the `sph_subs.h` library. As it occurs with the
 * single sphere case functions, in most cases, the results of calculations do
 * not fall back to fundamental data types. They are rather multi-component
 * structures. In order to manage access to such variety of return values, most
 * functions are declared as `void` and they operate on output arguments passed
 * by reference.
 */

#ifndef INCLUDE_TRA_SUBS_H_
#define INCLUDE_TRA_SUBS_H_
#endif

// Structures for TRAPPING
/*! \brief CIL data structure.
 *
 * A structure containing field expansion order configuration.
 */
struct CIL {
  //! Maximum external field expansion order
  int le;
  //! NLEM = LE * (LE + 2)
  int nlem;
  //! NLEMT = 2 * NLEM
  int nlemt;
  //! MXMPO is read from T-matrix with IS < 0 (not used here).
  int mxmpo;
  //! 2 * mxmpo - 1.
  int mxim;
};

/*! \brief CCR data structure.
 *
 * A structure containing geometrical asymmetry parameter normalization coefficients.
 */
struct CCR {
  //! First coefficient.
  double cof;
  //! Second coefficient.
  double cimu;
};
//End of TRAPPING structures

/*! C++ porting of CAMP
 *
 * \param ac: `complex double *` ANNOTATION: vector of complex amplitudes.
 * \param am0m: `complex double **` ANNOTATION: T-matrix of the particle.
 * \param ws: `complex double *` ANNOTATION: scattering amplitudes.
 * \param cil: `CIL *` Pointer to a CIL structure.
 *
 * This function builds the AC vector using AM0M and WS.
 */
void camp(dcomplex *ac, dcomplex **am0m, dcomplex *ws, CIL *cil);

/*! C++ porting of CZAMP
 *
 * \param ac: Vector of complex. ANNOTATION: vector of complex amplitudes.
 * \param amd: Matrix of complex. ANNOTATION: T-matrix.
 * \param indam: `int **`. ANNOTATION: map of matrix indexes.
 * \param ws: Vector of complex. ANNOTATION: scattering amplitudes.
 * \param cil: `CIL *` Pointer to a CIL structure.
 *
 * This function builds the AC vector using AMD, INDAM and WS.
 */
void czamp(dcomplex *ac, dcomplex **amd, int **indam, dcomplex *ws, CIL *cil);

/*! C++ porting of FFRF
 *
 * \param zpv: `double ****`. ANNOTATION: asymmetry tensor.
 * \param ac: `complex double *` ANNOTATION: vector of complex amplitudes.
 * \param ws: `complex double *` ANNOTATION: scattering amplitudes.
 * \param fffe: `double *`. ANNOTATION: extinction contributions to force.
 * \param fffs: `double *`. ANNOTATION: scattering contributions to force.
 * \param cil: `CIL *` Pointer to a CIL structure.
 * \param ccr: `CCR *` Pointer to a CCR structure.
 */
void ffrf(
	  double ****zpv, dcomplex *ac, dcomplex *ws, double *fffe,
	  double *fffs, CIL *cil, CCR *ccr
);

/*! C++ porting of FFRT
 *
 * \param ac: `complex double *` ANNOTATION: vector of complex amplitudes.
 * \param ws: `complex double *` ANNOTATION: scattering amplitudes.
 * \param ffte: `double *`. ANNOTATION: extinction contributions to torque.
 * \param ffts: `double *`. ANNOTATION: scattering contributions to torque.
 * \param cil: `CIL *` Pointer to a CIL structure.
 */
void ffrt(dcomplex *ac, dcomplex *ws, double *ffte, double *ffts, CIL *cil);

/*! C++ porting of FRFMER
 *
 * \param nkv: `int` ANNOTATION: number of wave vectors.
 * \param vkm: `double` ANNOTATION: wave number.
 * \param vknmx: `double` ANNOTATION: maximum wave number.
 * \param apfafa: `double` ANNOTATION: normalized aperture.
 * \param tra: `double` ANNOTATION: lens transmission.
 * \param spd: `double` ANNOTATION: cover thickness.
 * \param rir: `double` ANNOTATION: external over oil refractive index ratio.
 * \param ftcn: `double` ANNOTATION: refraction factor.
 * \param le: `int` ANNOTATION: external multipole truncation order.
 * \param lmode: `int` ANNOTATION: laser mode.
 * \param pmf: `double` ANNOTATION: aperture correction.
 * \param tt1: `Swap1 *` Pointer to first swap object.
 * \param tt2: `Swap2 *` Pointer to second swap object.
 * \return wk: `complex double *` ANNOTATION: incident amplitudes.
 */
dcomplex *frfmer(
	    int nkv, double vkm, double vknmx, double apfafa, double tra,
	    double spd, double rir, double ftcn, int le, int lmode, double pmf,
	    Swap1 *tt1, Swap2 *tt2
);

/*! C++ porting of PWMALP
 *
 * \param w: `complex double *` ANNOTATION: amplitudes.
 * \param up: `double *` ANNOTATION: unit polarization vector.
 * \param un: `double *` ANNOTATION: unit normal vector.
 * \param ylm: `complex double *` Field vector spherical harmonics.
 * \param lw: `int` ANNOTATION: order of amplitudes.
 */
void pwmalp(dcomplex **w, double *up, double *un, dcomplex *ylm, int lw);

/*! C++ porting of SAMP
 *
 * \param ac: `complex double *` ANNOTATION: vector of complex amplitudes.
 * \param tmsm: `complex double *` ANNOTATION: T-matrix M components.
 * \param tmse: `complex double *` ANNOTATION: T-matrix E components.
 * \param ws: `complex double *` ANNOTATION: scattering amplitudes.
 * \param cil: `CIL *` Pointer to a CIL structure.
 *
 * This function builds the AC vector using TMSM, TMSE and WS.
 */
void samp(dcomplex *ac, dcomplex *tmsm, dcomplex *tmse, dcomplex *ws, CIL *cil);

/*! C++ porting of SAMPOA
 *
 * \param ac: `complex double *` ANNOTATION: vector of complex amplitudes.
 * \param tms: `complex double **` ANNOTATION: T-matrix.
 * \param ws: `complex double *` ANNOTATION: scattering amplitudes.
 * \param cil: `CIL *` Pointer to a CIL structure.
 *
 * This function builds the AC vector using TMS and WS.
 */
void sampoa(dcomplex *ac, dcomplex **tms, dcomplex *ws, CIL *cil);

/*! C++ porting of WAMFF
 *
 * \param wk: `complex double *` ANNOTATION: incident amplitudes.
 * \param x: `double` ANNOTATION: X coordinate.
 * \param y: `double` ANNOTATION: Y coordinate.
 * \param z: `double` ANNOTATION: Z coordinate.
 * \param lm: `int` ANNOTATION: maximum multipole truncation order.
 * \param apfafa: `double` ANNOTATION: Normalized aperture.
 * \param tra: `double` ANNOTATION: Lens transmission.
 * \param spd: `double` ANNOTATION: Cover slip thickness.
 * \param rir: `double` ANNOTATION: External over oil refractive index ratio.
 * \param ftcn: `double` ANNOTATION: Refraction factor.
 * \param lmode: `int` ANNOTATION: Laser mode.
 * \param pmf: `double` ANNOTATION: aperture correction.
 */
void wamff(
	   dcomplex *wk, double x, double y, double z, int lm, double apfafa,
	   double tra, double spd, double rir, double ftcn, int lmode, double pmf
);
