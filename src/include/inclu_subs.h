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

/*! \file inclu_subs.h
 *
 * \brief C++ porting of INCLU functions and subroutines.
 *
 * This library includes a collection of functions that are used to solve the
 * scattering problem in the case of a sphere with a cluster of inclusions. Like
 * the other use cases, many of the functions declared here execute various
 * calculations on different data structures. In order to manage access to such
 * variety of calculations, most functions are declared as `void` and they operate
 * on output arguments passed by reference.
 */

#ifndef INCLUDE_INCLU_SUBS_H_
#define INCLUDE_INCLU_SUBS_H_

/*! \brief C++ porting of CNF.
 *
 * \param n: `int` Bessel y function order.
 * \param z: `dcomplex` Argument of Bessel y function.
 * \param nm: `int` Maximum computed order.
 * \param csj: `dcomplex *` ANNOTATION: Complex spherical J vector.
 * \param csy: `dcomplex *` Complex spherical Bessel functions up to desired order.
 */
void cnf(int n, dcomplex z, int nm, dcomplex *csj, dcomplex *csy);

/*! \brief C++ porting of EXMA.
 *
 * \param am: `dcomplex **` Field transition coefficients matrix.
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void exma(dcomplex **am, ParticleDescriptor *c1);

/*! \brief C++ porting of INCMS.
 *
 * \param am: `dcomplex **` Field transition coefficients matrix.
 * \param enti: `double` ANNOTATION: imaginary part of coating dielectric function.
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void incms(dcomplex **am, double enti, ParticleDescriptor *c1);

/*! \brief C++ porting of INDME.
 *
 * \param i: `int` 1-based sphere configuration index.
 * \param npnt: `int` ANNOTATION: Number of non transition layer integration points.
 * \param npntts: `int` ANNOTATION: Number of transition layer integrtion points.
 * \param vk: `double` Vacuum wave vector magnitude.
 * \param ent: `dcomplex` ANNOTATION: coating dielectric function.
 * \param enti: `double` ANNOTATION: imaginary part of coating dielectric function.
 * \param entn: `dcomplex` ANNOTATION: coating refractive index.
 * \param jer: `int &` Error code flag.
 * \param lcalc: `int &` Maximum order achieved in calculation.
 * \param arg: `dcomplex &` ANNOTATION: argument of calculation.
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void indme(
	   int i, int npnt, int npntts, double vk, dcomplex ent, double enti,
	   dcomplex entn, int &jer, int &lcalc, dcomplex &arg, ParticleDescriptor *c1
);

/*! \brief C++ porting of INSTR.
 *
 * \param rcf: `double **` Pointer to matrix of fractional radii.
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 */
void instr(double **rcf, ParticleDescriptor *c1);

/*! \brief C++ porting of OSPV.
 *
 * \param c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 * \param vk: `double` ANNOTATION: vacuum wave number.
 * \param sze: `double` ANNOTATION: size.
 * \param exri: `double` External medium refractive index.
 * \param entn: `dcomplex` Outer sphere refractive index.
 * \param enti: `double` Imaginary part of the outer medium refractive index.
 * \param jer: `int &` Reference to an integer error flag.
 * \param lcalc: `int &` Maximum order achieved in calculation.
 * \param arg: `dcomplex` Complex Bessel function argument.
 */
void ospv(ParticleDescriptor *c1, double vk, double sze, double exri, dcomplex entn, double enti, int &jer, int &lcalc, dcomplex &arg);

#endif // INCLUDE_INCLU_SUBS_H_
