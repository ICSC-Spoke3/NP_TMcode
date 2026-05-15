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

/**
 * \file magma_calls.h
 *
 * \brief C++ interface to MAGMA calls.
 *
 */

#ifndef INCLUDE_MAGMA_CALLS_H_
#define INCLUDE_MAGMA_CALLS_H_

#define MAX_GPU_STACK_DOUBLE 128

#ifdef USE_TARGET_OFFLOAD
#pragma omp begin declare target device_type(any)
/**
 * \brief C++ porting of CGEV. Get weight of T-matrix element (MAGMA GPU version).
 *
 * \param[in] ipamo: `int`
 * \param[in] mu: `int`
 * \param[in] l: `int`
 * \param[in] m: `int`
 * \return result: `double`
 */
double magma_cgev(int ipamo, int mu, int l, int m);
#pragma omp end declare target

/**
 * \brief Initialize a MAGMA compliant T-matrix leaving it on the GPU.
 *
 * \param[in,out] vec_am: `magmaDoubleComplex *` Vector form of the matrix.
 * \param[in] c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 * \param[in] rxx: `double *` Vector of sphere X coordinates.
 * \param[in] ryy: `double *` Vector of sphere Y coordinates.
 * \param[in] rzz: `double *` Vector of sphere Z coordinates.
 * \param[in] ind3j: `int *` Vector form of the 3J index look-up table.
 * \param[in] v3j0: `double *`
 * \param[in] vh: `magmaDoubleComplex *`
 * \param[in] vyhj: `magmaDoubleComplex *`
 * \param[in] vj0: `magmaDoubleComplex *`
 * \param[in] vyj0: `magmaDoubleComplex *`
 * \param[in] rmi: `magmaDoubleComplex *` Vector of Mie magnetic coefficients.
 * \param[in] rei: `magmaDoubleComplex *` Vector of Mie electric coefficients.
 * \param[in] device_id: `int` ID of the device to build the matrix on.
 * \return result: `magma_int_t` A return code (MAGMA_SUCCESS, if everything is fine).
 */
magma_int_t magma_cms(
  magmaDoubleComplex *vec_am, ParticleDescriptor *c1, double *rxx, double *ryy, double *rzz,
  int *ind3j, double *v3j0, magmaDoubleComplex *vh, magmaDoubleComplex *vyhj,
  magmaDoubleComplex *vj0, magmaDoubleComplex *vyj0, magmaDoubleComplex *rmi,
  magmaDoubleComplex* rei, int device_id
);

/**
 * \brief Compute the transfer vector from N2 to N1.
 *
 * This function computes the transfer vector going from N2 to N1, using either
 * Hankel, Bessel or Bessel from origin functions.
 *
 * \tparam IHI: `int` Function mode.
 * \param[in] ipamo: `int`
 * \param[in] nbl: `int` Block identifier.
 * \param[in] l1: `int` First L quantum number.
 * \param[in] m1: `int` First M quantum number.
 * \param[in] l2: `int` Second L quantum number.
 * \param[in] m2: `int` Second M quantum number.
 * \param[in] rxx: `double *` X coordinates of the spheres.
 * \param[in] ryy: `double *` Y coordinates of the spheres.
 * \param[in] rzz: `double *` Z coordinates of the spheres.
 * \param[in] ind3j: `int *` J vector index look-up table.
 * \param[in] v3j0: `double *` J0 vectors.
 * \param[in] vh: `magmaDoubleComplex *`
 * \param[in] vyhj: `magmaDoubleComplex *`
 * \param[in] vj0: `magmaDoubleComplex *`
 * \param[in] vyj0: `magmaDoubleComplex *`
 * \param[in] vj: `magmaDoubleComplex`
 * \param[in] li: `int` Maximum internal expansion order.
 * \param[in] le: `int` Maximum external expansion order.
 * \param[in,out] rac3j: `double[]` J connection vector.
 * \return result: `magmaDoubleComplex` Matrix element.
 */
#pragma omp begin declare target device_type(any)
template<int IHI> magmaDoubleComplex magma_ghit(
  int ipamo, int nbl, int l1, int m1, int l2, int m2,
  double *rxx, double *ryy, double *rzz, int *ind3j,
  double *v3j0, magmaDoubleComplex *vh, magmaDoubleComplex *vyhj,
  magmaDoubleComplex *vj0, magmaDoubleComplex *vyj0,
  magmaDoubleComplex vj, int li, int le, double rac3j[]
);
#pragma omp end declare target

#pragma omp begin declare target device_type(any)
/**
 * \brief Compute the 3j symbol for Clebsch-Gordan coefficients for JJ transitions.
 *
 * This function calculates the 3j(J,J2,J3;-M2-M3,M2,M3) symbol for the Clebsch-Gordan
 * coefficients. See Appendix a.3.1 in Borghese, Denti & Saija (2007). This function
 * is a specialized version to work on GPUs with MAGMA.
 *
 * \param[in] j2: `int` Value of J2.
 * \param[in] j3: `int` Value of J3.
 * \param[in] m2: `int` Value of M2.
 * \param[in] m3: `int` Value of M3.
 * \param[out] rac3j: `double[]` Vector of 3j symbols.
 */
void magma_r3jjr(int j2, int j3, int m2, int m3, double rac3j[]);
#pragma omp end declare target

/**
 * \brief Invert a pre-existing GPU complex matrix with double precision elements.
 *
 * \param[in,out] mat: Matrix of complex. The matrix to be inverted.
 * \param[in] n: `magma_int_t` The number of rows and columns of the [n x n] matrix.
 * \param[out] jer: `int &` Reference to an integer return flag.
 * \param[in] queue: `magma_queue_t` The MAGMA device communication queue.
 * \param[in] device_id: `int` ID of the device for matrix inversion offloading (optional, default 0).
 * \param[in] rs: `const RuntimeSettings &` Runtime settings instance (optional).
 */
void magma_zinvert_resident(
  magmaDoubleComplex *mat, magma_int_t n, int &jer, magma_queue_t queue, int device_id=0,
  const RuntimeSettings& rs=RuntimeSettings()
);

/**
 * \brief Compute the monocentered T-matrix
 *
 * \param[in] vec_am: `magmaDoubleComplex *` Vector form of the multi-centered matrix.
 * \param[in] c1: `ParticleDescriptor *` Pointer to a ParticleDescriptor instance.
 * \param[in] rxx: `double *` Vector of sphere X coordinates.
 * \param[in] ryy: `double *` Vector of sphere Y coordinates.
 * \param[in] rzz: `double *` Vector of sphere Z coordinates.
 * \param[in] ind3j: `int *` Vector form of the 3J index look-up table.
 * \param[in] v3j0: `double *`
 * \param[in] vh: `magmaDoubleComplex *`
 * \param[in] vyhj: `magmaDoubleComplex *`
 * \param[in] vj0: `magmaDoubleComplex *`
 * \param[in] vyj0: `magmaDoubleComplex *`
 * \param[in] sam_v: `magmaDoubleComplex *` Vector form of sums from vec_am.
 * \param[in] gis_v: `magmaDoubleComplex *`
 * \param[in] gls_v: `magmaDoubleComplex *`
 * \param[in] device_id: `int` ID of the device to build the matrix on.
 * \return result: `magma_int_t` A return code (MAGMA_SUCCESS, if everything is fine).
 */
magma_int_t magma_ztm(
  magmaDoubleComplex *vec_am, ParticleDescriptor *c1, double *rxx, double *ryy,
  double *rzz, int *ind3j, double *v3j0, magmaDoubleComplex *vh, magmaDoubleComplex *vyhj,
  magmaDoubleComplex *vj0, magmaDoubleComplex *vyj0, magmaDoubleComplex *sam_v,
  magmaDoubleComplex *gis_v, magmaDoubleComplex *gls_v, magmaDoubleComplex *vec_am0m,
  int device_id
);
#endif // USE_TARGET_OFFLOAD

/**
 * \brief Invert a complex matrix with double precision elements.
 *
 * This function is a first implementation of GPU use. The data transfer
 * steps are handled internally, via standard MAGMA calls, so that the matrix
 * to be inverted is copied from the host memory to the GPU, inverted, possibly
 * refined and finally taken back to the host, leaving the GPU free.
 *
 * \param[in,out] mat: Matrix of complex. The matrix to be inverted.
 * \param[in] n: `np_int` The number of rows and columns of the [n x n] matrix.
 * \param[out] jer: `int &` Reference to an integer return flag.
 * \param[in] device_id: `int` ID of the device for matrix inversion offloading (optional, default 0).
 * \param[in] rs: `const RuntimeSettings &` Runtime settings instance (optional).
 */
void magma_zinvert(
  dcomplex **mat, np_int n, int &jer, int device_id=0,
  const RuntimeSettings& rs=RuntimeSettings()
);

/**
 * \brief Perform Newton-Schulz iterative refinement of matrix inversion.
 *
 * In this function the residual of the inversion of a matrix A is evaluated as:
 *
 * R = A^-1 A - I
 *
 * and the convergence of refinement is estimated through the largest element
 * modulus left in R. Since Newton-Schultz iterative refinement is dangerous for
 * unstable matrices, the function first creates a host copy of the unrefined
 * inverted matrix and then attempts refinement. If the residue generated by
 * refinement does not improve over the original residue, the function returns
 * an error code that informs the calling function about corruption danger and
 * instructs it to use the back-up. In case of successful refinement, instead,
 * the calling function will be responsible of getting the refined matrix from
 * device back to the host.
 *
 * \param[in] rs: `const RuntimeSettings &` Runtime settings instance.
 * \param[in] a: `magmaDoubleComplex *` Pointer to the first element of the non-inverted matrix on host.
 * \param[in] m: `const magma_int_t` Number of rows / columns in a.
 * \param[in,out] d_a: `magmaDoubleComplex *` Pointer to the matrix on the GPU.
 * \param[in] queue: `magma_queue_t` GPU communication queue.
 * \return err: `magma_int_t` An error code (MAGMA_SUCCESS, if everything was fine).
 */
magma_int_t magma_newton(
  const RuntimeSettings& rs, magmaDoubleComplex* a, const magma_int_t m,
  magmaDoubleComplex* d_a, magma_queue_t queue
);

#endif
