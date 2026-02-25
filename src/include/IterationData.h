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

/*! \file IterationData.h
 *
 * \brief Multi-process communication data structures.
 *
 */

#ifndef INCLUDE_ITERATION_DATA_H_
#define INCLUDE_ITERATION_DATA_H_

// >>> DEFINITION OF ClusterIterationData CLASS <<<
/*! \brief A data structure representing the information used for a single scale
 * of the CLUSTER case.
 */
class ClusterIterationData {
public:
  //! \brief Pointer to a ParticleDescriptor structure.
  ParticleDescriptor *c1;
  //! \brief Vector of geometric asymmetry factors.
  double *gaps;
  //! \brief Components of extinction contribution to radiation torque on a single sphere along k.
  double **tqse;
  //! \brief Components of polarized extinction contribution to radiation torque on a single sphere along k.
  dcomplex **tqspe;
  //! \brief Components of scattering contribution to radiation torque on a single sphere along k.
  double **tqss;
  //! \brief Components of polarized scattering contribution to radiation torque on a single sphere along k.
  dcomplex **tqsps;
  //! \brief L-dependent coefficients of the geometric asymmetry parameter.
  double ****zpv;
  //! \brief Mean geometric asymmetry parameters.
  double **gapm;
  //! \brief Mean geometric asymmetry parameters referred to polarization plane.
  dcomplex **gappm;
  //! \brief Imaginary part of the harmonic functions argument.
  double *argi;
  //! \brief Argument of the harmonic functions referred to the scattering plane.
  double *args;
  //! \brief Geometric asymmetry parameters.
  double **gap;
  //! \brief Geometric asymmetry parameters referred to polarization plane.
  dcomplex **gapp;
  //! \brief Components of extinction contribution to radiation torque on the cluster along k.
  double **tqce;
  //! \brief Components of extinction contribution to radiation torque on the cluster along k referred to polarization plane.
  dcomplex **tqcpe;
  //! \brief Components of scattering contribution to radiation torque on the cluster along k.
  double **tqcs;
  //! \brief Components of scattering contribution to radiation torque on the cluster along k referred to polarization plane.
  dcomplex **tqcps;
  //! \brief Variation of unitary radiation vector.
  double *duk;
  //! \brief Cluster extinction cross-section components referred to scattering plane.
  double **cextlr;
  //! \brief Cluster extinction cross-section components referred to meridional plane.
  double **cext;
  //! \brief Cluster Mueller Transformation Matrix components referred to scattering plane.
  double **cmullr;
  //! \brief Cluster Mueller Transformation Matrix components referred to meridional plane.
  double **cmul;
  //! \brief Geometric asymmetry parameter components.
  double *gapv;
  //! \brief Radiation extinction torque components.
  double *tqev;
  //! \brief Radiation scattering torque components.
  double *tqsv;
  //! \brief Incident unitary vector components.
  double *u;
  //! \brief Scattered unitary vector components.
  double *us;
  //! \brief Normal unitary vector components.
  double *un;
  //! \brief Normal scattered unitary vector components.
  double *uns;
  //! \brief Incident unitary vector components on polarization plane.
  double *up;
  //! \brief Scattered unitary vector components on polarization plane.
  double *ups;
  //! \brief Mean unitary vector components normal to polarization plane.
  double *unmp;
  //! \brief Mean scattered unitary vector components normal to polarization plane.
  double *unsmp;
  //! \brief Mean incident unitary vector components on polarization plane.
  double *upmp;
  //! \brief Mean scattered unitary vector components on polarization plane.
  double *upsmp;
  //! \brief Scattering angle.
  double scan;
  //! \brief Control parameter on incidence direction referred to meridional plane.
  double cfmp;
  //! \brief Control parameter on scattering direction referred to meridional plane.
  double sfmp;
  //! \brief Control parameter on incidence direction referred to scattering plane.
  double cfsp;
  //! \brief Control parameter on scattering direction referred to scattering plane.
  double sfsp;
  //! \brief SQSFI = XI^-2
  double sqsfi;
  //! \brief Vectorized scattering coefficient matrix.
  dcomplex *am_vector;
  //! \brief Scattering coefficient matrix.
  dcomplex **am;
  //! \brief ANNOTATION: Argument of harmonic functions.
  dcomplex arg;
  //! \brief Vacuum magnitude of wave vector.
  double vk;
  //! \brief Wave number.
  double wn;
  //! \brief ANNOTATION: Normalization scale.
  double xip;
  //! \brief Number of scales (wavelengths) to be computed.
  int number_of_scales;
  //! \brief Size of the block of scales handled by the current process.
  int xiblock;
  //! \brief Index of the first scale handled by the current process.
  int firstxi;
  //! \brief Index of the last scale handled by the current process.
  int lastxi;
  //! \brief ID of the GPU used by one MPI process.
  int proc_device;
  //! \brief Refinement mode selction flag.
  int refinemode;
  //! \brief flag defining a first iteration
  bool is_first_scale;
  //! \brief Maximum number of refinement iterations.
  int maxrefiters;
  //! \brief Required accuracy level.
  double accuracygoal;

  /*! \brief `ClusterIterationData` default instance constructor.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \param mpidata: `mixMPI *` Pointer to a `mixMPI` object.
   * \param device_count: `const int` Number of offload devices available on the system.
   */
  ClusterIterationData(GeometryConfiguration *gconf, ScattererConfiguration *sconf, const mixMPI *mpidata, const int device_count);
  
  /*! \brief `ClusterIterationData` copy constructor.
   *
   * \param rhs: `const ClusterIterationData &` Reference to the `ClusterIterationData` object to be copied.
   */
  ClusterIterationData(const ClusterIterationData& rhs);

#ifdef MPI_VERSION
  /*! \brief `ClusterIterationData` MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to a `mixMPI` instance.
   * \param device_count: `const int` Number of offload devices available on the system.
   */
  ClusterIterationData(const mixMPI *mpidata, const int device_count);

  /*! \brief Broadcast over MPI the ClusterIterationData instance from MPI process 0 to all others.
   *
   * When using MPI, the initial ClusterIterationData instance created by MPI process 0
   * needs to be replicated on all other processes. This function sends it using
   * MPI broadcast calls. The MPI broadcast calls in this function must match those
   * in the constructor using the mixMPI pointer.
   *
   * \param mpidata: `mixMPI *` Pointer to the mpi structure used to do the MPI broadcast.
   */
  void mpibcast(const mixMPI *mpidata);
#endif // MPI_VERSION

  /*! \brief `ClusterIterationData` instance destroyer.
   */
  ~ClusterIterationData();

  /*! \brief Compute the memory requirements of an instance.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \return result: `long` Estimated size in bytes.
   */
  static long get_size(GeometryConfiguration *gconf, ScattererConfiguration *sconf);

  /*! \brief Update field expansion orders.
   *
   * \param rcf: `double **` Matrix of sphere fractional radii.
   * \param inner_order: `int` The new inner expansion order to be set.
   * \param outer_order: `int` The new outer expansion order to be set.
   * \return result: `int` An exit code (0 if successful).
   */
  int update_orders(double** rcf, int inner_order, int outer_order);
};
// >>> END OF ClusterIterationData CLASS DEFINITION <<<

// >>> DEFINITION OF InclusionIterationData CLASS <<<
/*! \brief A data structure representing the information used for a single scale
 * of the INCLUSION case.
 */
class InclusionIterationData {
protected:
  //! \brief Vectorized geometric asymmetry parameter components.
  double *vec_zpv;
  
public:
  //! \brief External layer index.
  int nimd;
  //! \brief External layer radius.
  double extr;
  
  //! \brief Pointer to a ParticleDescriptor structure.
  ParticleDescriptor *c1;
  //! \brief Vector of geometric asymmetry factors.
  double *gaps;
    //! \brief Components of extinction contribution to radiation torque on a single sphere along k.
  double **tqse;
  //! \brief Components of polarized extinction contribution to radiation torque on a single sphere along k.
  dcomplex **tqspe;
  //! \brief Components of scattering contribution to radiation torque on a single sphere along k.
  double **tqss;
  //! \brief Components of polarized scattering contribution to radiation torque on a single sphere along k.
  dcomplex **tqsps;
  //! \brief L-dependent coefficients of the geometric asymmetry parameter.
  double ****zpv;
  //! \brief Mean geometric asymmetry parameters.
  double **gapm;
  //! \brief Mean geometric asymmetry parameters referred to polarization plane.
  dcomplex **gappm;
  //! \brief Imaginary part of the harmonic functions argument.
  double *argi;
  //! \brief Argument of the harmonic functions referred to the scattering plane.
  double *args;
  //! \brief Geometric asymmetry parameters.
  double **gap;
  //! \brief Geometric asymmetry parameters referred to polarization plane.
  dcomplex **gapp;
  //! \brief Components of extinction contribution to radiation torque on the cluster along k.
  double **tqce;
  //! \brief Components of extinction contribution to radiation torque on the cluster along k referred to polarization plane.
  dcomplex **tqcpe;
  //! \brief Components of scattering contribution to radiation torque on the cluster along k.
  double **tqcs;
  //! \brief Components of scattering contribution to radiation torque on the cluster along k referred to polarization plane.
  dcomplex **tqcps;
  //! \brief Variation of unitary radiation vector.
  double *duk;
  //! \brief Cluster extinction cross-section components referred to scattering plane.
  double **cextlr;
  //! \brief Cluster extinction cross-section components referred to meridional plane.
  double **cext;
  //! \brief Cluster Mueller Transformation Matrix components referred to scattering plane.
  double **cmullr;
  //! \brief Cluster Mueller Transformation Matrix components referred to meridional plane.
  double **cmul;
  //! \brief Geometric asymmetry parameter components.
  double *gapv;
  //! \brief Radiation extinction torque components.
  double *tqev;
  //! \brief Radiation scattering torque components.
  double *tqsv;
  //! \brief Incident unitary vector components.
  double *u;
  //! \brief Scattered unitary vector components.
  double *us;
  //! \brief Normal unitary vector components.
  double *un;
  //! \brief Normal scattered unitary vector components.
  double *uns;
  //! \brief Incident unitary vector components on polarization plane.
  double *up;
  //! \brief Scattered unitary vector components on polarization plane.
  double *ups;
  //! \brief Mean unitary vector components normal to polarization plane.
  double *unmp;
  //! \brief Mean scattered unitary vector components normal to polarization plane.
  double *unsmp;
  //! \brief Mean incident unitary vector components on polarization plane.
  double *upmp;
  //! \brief Mean scattered unitary vector components on polarization plane.
  double *upsmp;
  //! \brief Scattering angle.
  double scan;
  //! \brief Control parameter on incidence direction referred to meridional plane.
  double cfmp;
  //! \brief Control parameter on scattering direction referred to meridional plane.
  double sfmp;
  //! \brief Control parameter on incidence direction referred to scattering plane.
  double cfsp;
  //! \brief Control parameter on scattering direction referred to scattering plane.
  double sfsp;
  //! \brief SQSFI = XI^-2
  double sqsfi;
  //! \brief Vectorized scattering coefficient matrix.
  dcomplex *am_vector;
  //! \brief Scattering coefficient matrix.
  dcomplex **am;
  //! \brief Argument of harmonic functions.
  dcomplex arg;
  //! \brief Vacuum magnitude of wave vector.
  double vk;
  //! \brief Wave number.
  double wn;
  //! \brief ANNOTATION: Normalization scale.
  double xip;
  //! \brief Number of scales (wavelengths) to be computed.
  int number_of_scales;
  //! \brief Size of the block of scales handled by the current process.
  int xiblock;
  //! \brief Index of the first scale handled by the current process.
  int firstxi;
  //! \brief Index of the last scale handled by the current process.
  int lastxi;
  //! \brief ID of the GPU used by one MPI process.
  int proc_device;
  //! \brief Refinement mode selction flag.
  int refinemode;
  //! \brief flag defining a first iteration
  bool is_first_scale;
  //! \brief Maximum number of refinement iterations.
  int maxrefiters;
  //! \brief Required accuracy level.
  double accuracygoal;

  /*! \brief `InclusionIterationData` default instance constructor.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \param mpidata: `mixMPI *` Pointer to a `mixMPI` object.
   * \param device_count: `const int` Number of offload devices available on the system.
   */
  InclusionIterationData(GeometryConfiguration *gconf, ScattererConfiguration *sconf, const mixMPI *mpidata, const int device_count);
  
  /*! \brief `InclusionIterationData` copy constructor.
   *
   * \param rhs: `const InclusionIterationData &` Reference to the `InclusionIterationData` object to be copied.
   */
  InclusionIterationData(const InclusionIterationData& rhs);

#ifdef MPI_VERSION
  /*! \brief `InclusionIterationData` MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to a `mixMPI` instance.
   * \param device_count: `const int` Number of offload devices available on the system.
   */
  InclusionIterationData(const mixMPI *mpidata, const int device_count);

  /*! \brief Broadcast over MPI the InclusionIterationData instance from MPI process 0 to all others.
   *
   * When using MPI, the initial InclusionIterationData instance created by MPI process 0
   * needs to be replicated on all other processes. This function sends it using
   * MPI broadcast calls. The MPI broadcast calls in this function must match those
   * in the constructor using the mixMPI pointer.
   *
   * \param mpidata: `mixMPI *` Pointer to the mpi structure used to do the MPI broadcast.
   */
  void mpibcast(const mixMPI *mpidata);
#endif

  /*! \brief `InclusionIterationData` instance destroyer.
   */
  ~InclusionIterationData();

  /*! \brief Compute the memory requirements of an instance.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \return result: `long` Estimated size in bytes.
   */
  static long get_size(GeometryConfiguration *gconf, ScattererConfiguration *sconf);

  /*! \brief Update field expansion orders.
   *
   * \param rcf: `double **` Matrix of sphere fractional radii.
   * \param inner_order: `int` The new inner expansion order to be set.
   * \param outer_order: `int` The new outer expansion order to be set.
   * \return result: `int` An exit code (0 if successful).
   */
  int update_orders(double** rcf, int inner_order, int outer_order);
};
// >>> END OF InclusionIterationData CLASS DEFINITION <<< //

// >>> DEFINITION OF SphereIterationData CLASS <<<
/*! \brief A data structure representing the information used for a single scale
 * of the SPHERE case.
 */
class SphereIterationData {
protected:
  //! \brief Number of spheres
  int _nsph;
  //! \brief Maximum field expansion order.
  int _lm;
  //! \brief Vector of Mueller matrix components.
  double *vec_cmul;
  //! \brief Vector of Mueller matrix components referred to meridional plane.
  double *vec_cmullr;
  //! Vectorized TQSPE.
  dcomplex *vec_tqspe;
  //! Vectorized TQSPS.
  dcomplex *vec_tqsps;
  //! Vectorized TQSE.
  double *vec_tqse;
  //! Vectorized TQSS.
  double *vec_tqss;
  //! Vectorized ZPV.
  double *vec_zpv;
  
public:
  //! \brief Cosine of incident radiation azimuth.
  double cost;
  //! \brief Sine of incident radiation azimuth.
  double sint;
  //! \brief Cosine of incident radiation elevation.
  double cosp;
  //! \brief Sine of incident radiation elevation.
  double sinp;
  //! \brief Cosine of scattered radiation azimuth.
  double costs;
  //! \brief Sine of scattered radiation azimuth.
  double sints;
  //! \brief Cosine of scattered radiation elevation.
  double cosps;
  //! \brief Sine of scattered radiation elevation.
  double sinps;
  //! \brief Vacuum magnitude of wave vector.
  double vk;
  //! \brief Wave number.
  double wn;
  //! \brief ANNOTATION: Normalization scale.
  double xip;
  //! \brief Number of scales (wavelengths) to be computed.
  int number_of_scales;
  //! \brief Size of the block of scales handled by the current process.
  int xiblock;
  //! \brief Index of the first scale handled by the current process.
  int firstxi;
  //! \brief Index of the last scale handled by the current process.
  int lastxi;
  //! \brief Argument of harmonic functions.
  dcomplex arg;
  //! \brief S0 = FSAS / (4 PI K^3).
  dcomplex s0;
  //! \brief Total forward scattering amplitude of the spheres.
  dcomplex tfsas;
  //! \brief Pointer to a sphere particle descriptor.
  ParticleDescriptor *c1;
  //! \brief Imaginary part of `arg`.
  double *argi;
  //! \brief `arg` squared.
  double *args;
  //! \brief Scattering angle.
  double scan;
  //! \brief Control parameter on incidence direction referred to meridional plane.
  double cfmp;
  //! \brief Control parameter on scattering direction referred to meridional plane.
  double sfmp;
  //! \brief Control parameter on incidence direction referred to scattering plane.
  double cfsp;
  //! \brief Control parameter on scattering direction referred to scattering plane.
  double sfsp;
  //! \brief Geometry asymmetry parameter for spheres.
  double *gaps;
  //! \brief Variation of unitary wave vector.
  double *duk;
  //! \brief Incidence direction unitary vector.
  double *u;
  //! \brief Scattering direction unitary vector.
  double *us;
  //! \brief Normal direction unitary vector.
  double *un;
  //! \brief Scattering normal direction unitary vector.
  double *uns;
  //! \brief Polarization direction unitary vector.
  double *up;
  //! \brief Scattered polarization direction unitary vector.
  double *ups;
  //! \brief Polarization direction unitary vector referred to meridional plane.
  double *upmp;
  //! \brief Scattered polarization direction unitary vector referred to meridional plane.
  double *upsmp;
  //! \brief Normal direction unitary vector referred to meridional plane.
  double *unmp;
  //! \brief Scattering normal direction unitary vector referred to meridional plane.
  double *unsmp;
  //! \brief Mueller matrix components.
  double **cmul;
  //! \brief Mueller matrix components referred to meridional plane.
  double **cmullr;
  //! \brief Polarization-dependent extinction contribution to torque for each sphere.
  dcomplex **tqspe;
  //! \brief Polarization-dependent scattering contribution to torque for each sphere.
  dcomplex **tqsps;
  //! \brief Extinction contribution to torque for each sphere.
  double **tqse;
  //! \brief Scattering contribution to torque for each sphere.
  double **tqss;
  //! \brief Scattering coefficients tensor.
  double ****zpv;
  //! \brief flag for first time initialisation
  bool is_first_scale;
  
  /*! \brief `SphereIterationData` default instance constructor.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \param mpidata: `mixMPI *` Pointer to a `mixMPI` object.
   * \param device_count: `const int` Number of offload devices available on the system.
   */
  SphereIterationData(GeometryConfiguration *gconf, ScattererConfiguration *sconf, const mixMPI *mpidata, const int device_count);
  
  /*! \brief `SphereIterationData` copy constructor.
   *
   * \param rhs: `const SphereIterationData &` Reference to the object to be copied.
   */
  SphereIterationData(const SphereIterationData& rhs);

#ifdef MPI_VERSION
  /*! \brief `SphereIterationData` MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to a `mixMPI` instance.
   * \param device_count: `const int` Number of offload devices available on the system.
   */
  SphereIterationData(const mixMPI *mpidata, const int device_count);

  /*! \brief Broadcast over MPI the `SphereIterationData` instance from MPI process 0 to all others.
   *
   * When using MPI, the initial InclusionIterationData instance created by
   * MPI process 0 needs to be replicated on all other processes. This
   * function sends it using MPI broadcast calls. The MPI broadcast calls in
   * this function must match those in the constructor using the mixMPI pointer.
   *
   * \param mpidata: `mixMPI *` Pointer to `mixMPI` instance.
   */
  int mpibcast(const mixMPI *mpidata);
#endif // MPI_VERSION

  /*! \brief `SphereIterationData` instance destroyer.
   */
  ~SphereIterationData();
  
  /*! \brief Update field expansion order.
   *
   * \param order: `int` The new expansion order to be set.
   * \return result: `int` An exit code (0 if successful).
   */
  int update_order(int order);
};
// >>> END OF SphereIterationData CLASS DEFINITION <<<

#endif // INCLUDE_ITERATION_DATA_H_
