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

/*! \file Commons.h
 *
 * \brief C++ porting of common data structures.
 *
 * Many functions of the original FORTRAN code share complex data blocks in
 * form of COMMON blocks. This poses the limit of freezing the structure of
 * the data blocks in the code, therefore implying the necessity to modify
 * the code to adapt it to the input and to recompile before running. C++,
 * on the contrary, offers the possibility to represent the necessary data
 * structures as classes that can dynamically instantiate the shared information
 * in the most convenient format for the current configuration. This approach
 * adds an abstraction layer that lifts the need to modify and recompile the
 * code depending on the input.
 *
 */

#ifndef INCLUDE_COMMONS_H_
#define INCLUDE_COMMONS_H_

class ParticleDescriptor;

/*! \brief Structure with essential MPI data.
 */
class mixMPI {
public:
  //! \brief was MPI initialised?
  bool mpirunning;
  //! \brief MPI rank
  int rank;
  //! \brief MPI nprocs
  int nprocs;

  /*! \brief empty mixMPI instance constructor.
   */
  mixMPI();

  /*! \brief mixMPI instance constructor from an actual MPI communicator.
   */
#ifdef MPI_VERSION
  mixMPI(MPI_Comm comm);
#endif
  
  /*! \brief mixMPI instance constructor copying its contents from a preexisting object.
   */
  mixMPI(const mixMPI& rhs);

  /*! \brief mixMPI instance destroyer.
   */
  ~mixMPI();
};

/*! \brief Basic data structure describing the particle model and its interaction with fields.
 *
 * This class forms a base of the data structure collections that are used by the
 * the code to handle different scattering problems. The class implements real
 * members that are used by all code sections.
 */
class ParticleDescriptor {
protected:
  // >>> COMMON TO ALL DESCRIPTOR TYPES <<< //
  //! \brief Sub-class identification code.
  short _class_type;
  //! \brief Number of spheres composing the model.
  int _nsph;
  //! \brief Maximum internal field expansion order.
  int _li;
  //! \brief Maximum number of layers in known sphere types.
  int _max_layers;
  //! \brief Number of different sphere types.
  int _num_configurations;
  //! \brief Total number of layers from all sphere types.
  int _num_layers;
  //! \brief NHSPO = 2 * MAX(NPNT,NPNTTS) - 1
  int _nhspo;
  //! \brief Number of points for numerical integration in layered spheres.
  int _npnt;
  //! \brief Number of points for numerical integration in transition layer.
  int _npntts;
  //! \brief Contiguous space for RMI.
  dcomplex *vec_rmi;
  //! \brief Contiguous space for REI.
  dcomplex *vec_rei;
  //! \brief Contiguous space for W.
  dcomplex *vec_w;
  //! \brief Contiguous space for RC.
  double *vec_rc;
  // >>> END OF SECTION COMMON TO ALL DESCRIPTOR TYPES <<< //

  // >>> NEEDED BY SPHERE AND CLUSTER <<< //
  //! \brief Contiguous space for SAS.
  dcomplex *vec_sas;
  //! \brief Contiguous space for VINTS.
  dcomplex *vec_vints;
  // >>> END OF SECTION NEEDED BY SPHERE AND CLUSTER <<< //

  // >>> NEEDED BY CLUSTER <<< //
  //! \brief Contiguous space for TSAS.
  dcomplex *vec_tsas;
  //! \brief Contiguous space for GIS.
  dcomplex *vec_gis;
  //! \brief Contiguous space for GLS.
  dcomplex *vec_gls;
  //! \brief Contiguous space for SAM.
  dcomplex *vec_sam;
  // >>> END OF SECTION NEEDED BY CLUSTER <<< //
  
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  //! \brief Maximum external field expansion order.
  int _le;
  //! \brief Maximum field expansion order.
  int _lm;
  //! \brief NLIM = LI * (LI + 2)
  int _nlim;
  //! \brief NLEM = LE * (LE + 2)
  int _nlem;
  //! \brief NLEMT = 2 * NLEM
  int _nlemt;
  //! \brief NCOU = NSPH * NSPH - 1
  int _ncou;
  //! \brief LITPO = 2 * LI + 1
  int _litpo;
  //! \brief LITPOS = LITPO * LITPO
  int _litpos;
  //! \brief LMPO = LM + 1
  int _lmpo;
  //! \brief LMTPO = LI + LE + 1
  int _lmtpo;
  //! \brief LMTPOS = LMTPO * LMTPO
  int _lmtpos;
  //! \brief NV3J = (LM * (LM + 1) * (2 * LM + 7)) / 6
  int _nv3j;
  //! \brief NDI = NSPH * NLIM
  int _ndi;
  //! \brief NDIT = 2 * NSPH * NLIM
  int _ndit;
  //! \brief Contiguous space for AM0M.
  dcomplex *vec_am0m;
  //! \brief Contiguous space for FSAC.
  dcomplex *vec_fsac;
  //! \brief Contiguous space for SAC.
  dcomplex *vec_sac;
  //! \brief Contiguous space for FSACM.
  dcomplex *vec_fsacm;
  //! \brief Contiguous space for IND3J.
  int *vec_ind3j;
  // >>> END OF SECTION NEEDED BY CLUSTER AND INCLU <<< //

  // >>> NEEDED BY INCLU <<< //
  //! \brief NDM = NDIT + NLEMT
  int _ndm;
  //! \brief Contiguous space for AT.
  dcomplex *vec_at;
  // >>> END OF SECTION NEEDED BY INCLU <<< //
public:
  // >>> COMMON TO ALL DESCRIPTOR TYPES <<< //
  //! \brief Base sub-class identification code.
  const static short BASE_TYPE = 0;
  //! \brief Sphere sub-class identification code.
  const static short SPHERE_TYPE = 1;
  //! \brief Cluster sub-class identification code.
  const static short CLUSTER_TYPE = 2;
  //! \brief Inclusion sub-class identification code.
  const static short INCLUSION_TYPE = 3;
  
  //! \brief Sub-class identification code.
  const short& class_type = _class_type;
  //! \brief Read-only view of number of spheres composing the model.
  const int &nsph = _nsph;
  //! \brief Read-only view of maximum internal field expansion order.
  const int &li = _li;
  //! \brief Read-only view of the maximum number of layers in types.
  const int &max_layers = _max_layers;
  //! \brief Read-only view of number of different sphere types.
  const int &num_configurations = _num_configurations;
  //! \brief Read-only view of total number of layers from all sphere types.
  const int &num_layers = _num_layers;
  //! \brief Read-only view of NHSPO.
  const int &nhspo = _nhspo;
  //! \brief Read-only view on number of points for numerical integration in layered spheres.
  const int &npnt = _npnt;
  //! \brief Read-only view on number of points for numerical integration in transition layer.
  const int &npntts = _npntts;
  //! \brief Matrix of Mie scattering b-coefficients for different orders and spheres.
  dcomplex **rmi;
  //! \brief Matrix of Matrix of Mie scattering a-coefficients for different orders and spheres.
  dcomplex **rei;
  //! \brief Matrix of multipole amplitudes for incident and scattered fields.
  dcomplex **w;
  //! \brief Vector of intensity components.
  dcomplex *vint;
  //! \brief Vector of sphere Cartesian X coordinates.
  double *rxx;
  //! \brief Vector of sphere Cartesian Y coordinates.
  double *ryy;
  //! \brief Vector of sphere Cartesian Z coordinates.
  double *rzz;
  //! \brief Vector of sphere radii.
  double *ros;
  //! \brief Matrix of sphere fractional transition radii.
  double **rc;
  //! \brief Vector of sphere type identifiers.
  int *iog;
  //! \brief Vector of number of layers in sphere type.
  int *nshl;
  //! \brief TBD
  dcomplex *ris;
  //! \brief TBD
  dcomplex *dlri;
  //! \brief Vector of dielectric constants.
  dcomplex *dc0;
  //! \brief Vector of complex refractive indices.
  dcomplex *vkt;
  //! \brief Vector of sizes in units of 2*PI/LAMBDA
  double *vsz;
  //! \brief Total geometric cross-section.
  double gcs;
  // >>> END OF SECTION COMMON TO ALL DESCRIPTOR TYPES <<< //
  
  // >>> NEEDED BY SPHERE AND CLUSTER <<< //
  //! \brief Tensor of sphere scattering amplitudes.
  dcomplex ***sas;
  //! \brief Array of vectors of intensity components for single spheres.
  dcomplex **vints;
  //! \brief Vector of sphere forward scattering amplitudes.
  dcomplex *fsas;
  //! \brief Vector of scattering cross-sections for spheres.
  double *sscs;
  //! \brief Vector of extinction cross-sections for spheres.
  double *sexs;
  //! \brief Vector of absorption cross-sections for spheres.
  double *sabs;
  //! \brief Vector of scattering efficiencies for spheres.
  double *sqscs;
  //! \brief Vector of extinction efficiencies for spheres.
  double *sqexs;
  //! \brief Vector of absorption efficiencies for spheres.
  double *sqabs;
  //! \brief Vector of geometric cross-sections for spheres.
  double *gcsv;
  // >>> END OF SECTION NEEDED BY SPHERE AND CLUSTER <<< //

  // >>> NEEDED BY CLUSTER <<<
  //! \brief Vector of field intensity components.
  dcomplex *vintt;
  //! \brief Total forward scattering amplitude of the spheres.
  dcomplex tfsas;
  //! \brief Total scattering amplitude of the spheres.
  dcomplex **tsas;
  //! \brief Sphere scattering cross-section.
  double scs;
  //! \brief Sphere extinction cross-section.
  double ecs;
  //! \brief Sphere absorption cross-section.
  double acs;
  //! \brief TBD.
  dcomplex **gis;
  //! \brief TBD.
  dcomplex **gls;
  //! \brief Mean scattering amplitude components.
  dcomplex **sam;
  // >>> END OF SECTION NEEDED BY CLUSTER <<< //
  
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  //! \brief Read-only view of maximum external field expansion order.
  const int& le = _le;
  //! \brief Read-only view of maximum field expansion order.
  const int& lm = _lm;
  //! \brief Read-only view of NLIM.
  const int& nlim = _nlim;
  //! \brief Read-only view of NLEM.
  const int& nlem = _nlem;
  //! \brief Read-only view of NLEMT.
  const int& nlemt = _nlemt;
  //! \brief Read-only view of NCOU.
  const int& ncou = _ncou;
  //! \brief Read-only view of LITPO.
  const int& litpo = _litpo;
  //! \brief Read-only view of LITPOS.
  const int& litpos = _litpos;
  //! \brief Read-only view of LMPO.
  const int& lmpo = _lmpo;
  //! \brief Read-only view of LMTPO.
  const int& lmtpo = _lmtpo;
  //! \brief Read-only view of LMTPOS.
  const int& lmtpos = _lmtpos;
  //! \brief Read-only view of NV3J.
  const int& nv3j = _nv3j;
  //! \brief Read-only view of NDI.
  const int& ndi = _ndi;
  //! \brief Read-only view of NDIT.
  const int& ndit = _ndit;
  //! \brief Read-only view of NDM.
  const int& ndm = _ndm;

  //! \brief TBD
  dcomplex *vh;
  //! \brief TBD
  dcomplex *vj0;
  //! \brief TBD
  dcomplex *vyhj;
  //! \brief TBD
  dcomplex *vyj0;
  //! \brief TBD
  dcomplex vj;
  //! \brief Transition matrix
  dcomplex **am0m;
  //! \brief Cluster forward scattering amplitude.
  dcomplex **fsac;
  //! \brief Cluster scattering amplitude.
  dcomplex **sac;
  //! \brief Mean cluster forward scattering amplitude.
  dcomplex **fsacm;
  //! \brief Mean intensity components vector.
  dcomplex *vintm;
  //! \brief Cluster polarized scattering cross-sections.
  dcomplex *scscp;
  //! \brief Cluster polarized extinction cross-sections.
  dcomplex *ecscp;
  //! \brief Mean cluster polarized scattering cross-sections.
  dcomplex *scscpm;
  //! \brief Mean cluster polarized extinction cross-sections.
  dcomplex *ecscpm;
  //! \brief TBD
  double *v3j0;
  //! \brief Cluster scattering cross-sections.
  double *scsc;
  //! \brief Cluster extinction cross-sections.
  double *ecsc;
  //! \brief Mean cluster scattering cross-sections.
  double *scscm;
  //! \brief Mean cluster extinction cross-sections.
  double *ecscm;
  //! \brief J-vector components index matrix.
  int **ind3j;
  //! \brief J-vector boundary values. QUESTION: correct?
  double *rac3j;
  // >>> END OF SECTION NEEDED BY CLUSTER AND INCLU <<< //

  // >>> NEEDED BY INCLU <<< //
  //! \brief TBD
  dcomplex *rm0;
  //! \brief TBD
  dcomplex *re0;
  //! \brief TBD
  dcomplex *rmw;
  //! \brief TBD
  dcomplex *rew;
  //! \brief TBD
  dcomplex *tm;
  //! \brief TBD
  dcomplex *te;
  //! \brief TBD
  dcomplex *tm0;
  //! \brief TBD
  dcomplex *te0;
  //! \brief TBD
  dcomplex **at;
  // >>> END OF SECTION NEEDED BY INCLU <<< //
  
  /*! \brief ParticleDescriptor instance constructor.
   *
   * \param gconf: `const GeometryConfiguration *` Pointer to GeometryConfiguration instance.
   * \param sconf: `const ScattererConfiguration *` Pointer to ScattererConfiguration instance.
   */
  ParticleDescriptor(GeometryConfiguration *gconf, ScattererConfiguration *sconf);

  /*! \brief ParticleDescriptor copy constructor.
   *
   * \param rhs: `const ParticleDescriptor &` Reference to ParticleDescriptor object to be copied.
   */
  ParticleDescriptor(const ParticleDescriptor &rhs);

#ifdef MPI_VERSION
  /*! \brief ParticleDescriptor MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  ParticleDescriptor(const mixMPI *mpidata);
  
  /*! \brief Broadcast ParticleDescriptor instance over MPI.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  void mpibcast(const mixMPI *mpidata);
#endif // MPI_VERSION
  
  /*! \brief ParticleDescriptor instance destroyer.
   */
  ~ParticleDescriptor();

  /*! \brief Interface function to return the descriptor type as string.
   *
   * \return descriptor_type: `string` The descriptor type name.
   */
  virtual std::string get_descriptor_type() { return "base descriptor"; }

  /*! \brief Compute the memory requirements of an instance.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \return result: `long` Estimated size in bytes.
   */
  static long get_size(GeometryConfiguration *gconf, ScattererConfiguration *sconf);
};

/*! \brief The data structure describing a particle model made by a cluster of spheres.
 *
 * This class is used to solve the problem of a cluster of spherical particles.
 * It replaces the C1 and C1_AddOns class implementation of the C1 FORTRAN common
 * block.
 */
class ParticleDescriptorCluster: public ParticleDescriptor {
public:
  /*! \brief ParticleDescriptorCluster instance constructor.
   *
   * \param gconf: `const GeometryConfiguration *` Pointer to GeometryConfiguration instance.
   * \param sconf: `const ScattererConfiguration *` Pointer to ScattererConfiguration instance.
   */
  ParticleDescriptorCluster(GeometryConfiguration *gconf, ScattererConfiguration *sconf);

  /*! \brief ParticleDescriptorCluster copy constructor.
   *
   * \param rhs: `const ParticleDescriptorCluster &` Reference to ParticleDescriptorCluster object to be copied.
   */
  ParticleDescriptorCluster(const ParticleDescriptorCluster &rhs);

#ifdef MPI_VERSION
  /*! \brief ParticleDescriptorCluster MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  ParticleDescriptorCluster(const mixMPI *mpidata);
  
  /*! \brief Broadcast ParticleDescriptorCluster instance over MPI.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  void mpibcast(const mixMPI *mpidata);
#endif // MPI_VERSION
  
  /*! \brief Interface function to return the descriptor type as string.
   *
   * \return descriptor_type: `string` The descriptor type name.
   */
  std::string get_descriptor_type() override { return "cluster descriptor"; }

  /*! \brief Compute the memory requirements of an instance.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \return result: `long` Estimated size in bytes.
   */
  static long get_size(GeometryConfiguration *gconf, ScattererConfiguration *sconf);
  
  /*! \brief Update the field expansion orders.
   *
   * \param inner_order: `int` The new inner expansion order to be set.
   * \param outer_order: `int` The new outer expansion order to be set.
   * \return result: `int` An exit code (0 if successful).
   */
  int update_orders(int inner_order, int outer_order);
};

/*! \brief The data structure describing a particle model for a sphere with inclusions.
 *
 * This class is used to solve the problem of a spherical particle with a
 * cluster of inclusions. It replaces the C1 and C1_AddOns class implementation
 * of the C1 FORTRAN common block.
 */
class ParticleDescriptorInclusion: public ParticleDescriptor {
public:
  /*! \brief ParticleDescriptorInclusion instance constructor.
   *
   * \param gconf: `const GeometryConfiguration *` Pointer to GeometryConfiguration instance.
   * \param sconf: `const ScattererConfiguration *` Pointer to ScattererConfiguration instance.
   */
  ParticleDescriptorInclusion(GeometryConfiguration *gconf, ScattererConfiguration *sconf);

  /*! \brief ParticleDescriptorInclusion copy constructor.
   *
   * \param rhs: `const ParticleDescriptorInclusion &` Reference to ParticleDescriptorInclusion object to be copied.
   */
  ParticleDescriptorInclusion(const ParticleDescriptorInclusion &rhs);

#ifdef MPI_VERSION
  /*! \brief ParticleDescriptorInclusion MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  ParticleDescriptorInclusion(const mixMPI *mpidata);
  
  /*! \brief Broadcast ParticleDescriptorInclusion instance over MPI.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  void mpibcast(const mixMPI *mpidata);
#endif // MPI_VERSION
  
  /*! \brief Interface function to return the descriptor type as string.
   *
   * \return descriptor_type: `string` The descriptor type name.
   */
  std::string get_descriptor_type() override { return "inclusion descriptor"; }

  /*! \brief Compute the memory requirements of an instance.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \return result: `long` Estimated size in bytes.
   */
  static long get_size(GeometryConfiguration *gconf, ScattererConfiguration *sconf);
  
  /*! \brief Update the field expansion orders.
   *
   * \param inner_order: `int` The new inner expansion order to be set.
   * \param outer_order: `int` The new outer expansion order to be set.
   * \return result: `int` An exit code (0 if successful).
   */
  int update_orders(int inner_order, int outer_order);
};

/*! \brief The data structure describing a spherical particle model.
 *
 * This class is used to solve the problem of a single spherical particle. It
 * replaces the old C1 class implementation of the corresponding FORTRAN common
 * block.
 */
class ParticleDescriptorSphere: public ParticleDescriptor {
public:
  /*! \brief ParticleDescriptorSphere instance constructor.
   *
   * \param gconf: `const GeometryConfiguration *` Pointer to GeometryConfiguration instance.
   * \param sconf: `const ScattererConfiguration *` Pointer to ScattererConfiguration instance.
   */
  ParticleDescriptorSphere(GeometryConfiguration *gconf, ScattererConfiguration *sconf);

  /*! \brief ParticleDescriptorSphere copy constructor.
   *
   * \param rhs: `const ParticleDescriptorSphere &` Reference to ParticleDescriptorSphere object to be copied.
   */
  ParticleDescriptorSphere(const ParticleDescriptorSphere &rhs);

#ifdef MPI_VERSION
  /*! \brief ParticleDescriptorSphere MPI constructor.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  ParticleDescriptorSphere(const mixMPI *mpidata);
  
  /*! \brief Broadcast ParticleDescriptorSphere instance over MPI.
   *
   * \param mpidata: `const mixMPI *` Pointer to mixMPI instance.
   */
  void mpibcast(const mixMPI *mpidata);
#endif // MPI_VERSION
  
  /*! \brief Interface function to return the descriptor type as string.
   *
   * \return descriptor_type: `string` The descriptor type name.
   */
  std::string get_descriptor_type() override { return "sphere descriptor"; }

  /*! \brief Compute the memory requirements of an instance.
   *
   * \param gconf: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` object.
   * \param sconf: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` object.
   * \return result: `long` Estimated size in bytes.
   */
  static long get_size(GeometryConfiguration *gconf, ScattererConfiguration *sconf);
  
  /*! \brief Update the field expansion order.
   *
   * \param order: `int` The new field expansion order to be set.
   * \return result: `int` An exit code (0 if successful).
   */
  int update_order(int order);
};

/*! \brief A data structure representing the angles to be evaluated in the problem.
 *
 */
class ScatteringAngles {
protected:
  //! \brief Number of incident field azimuth angles.
  int _nth;
  //! \brief Number of scattered field azimuth angles.
  int _nths;
  //! \brief Number of incident field elevation angles.
  int _nph;
  //! \brief Number of scattered field elevation angles.
  int _nphs;
  //! \brief Number of incident field propagation angles.
  int _nk;
  //! \brief Number of scattered field propagation angles.
  int _nks;
  //! \brief Total number of field propagation angles.
  int _nkks;
  //! \brief First incident field azimuth angle.
  double _th;
  //! \brief Incident field azimuth angle increment.
  double _thstp;
  //! \brief Last incident field azimuth angle.
  double _thlst;
  //! \brief First scattered field azimuth angle.
  double _ths;
  //! \brief Scattered field azimuth angle increment.
  double _thsstp;
  //! \brief Last scattered field azimuth angle.
  double _thslst;
  //! \brief First incident field elevation angle.
  double _ph;
  //! \brief Incident field elevation angle increment.
  double _phstp;
  //! \brief Last incident field elevation angle.
  double _phlst;
  //! \brief First scattered field elevation angle.
  double _phs;
  //! \brief Scattered field elevation angle increment.
  double _phsstp;
  //! \brief Last scattered field elevation angle.
  double _phslst;
  //! \brief Azimuth scattering deflection.
  double _thsca;

public:
  //! \brief Read only view of `_nth`.
  const int& nth = _nth;
  //! \brief Read only view of `_nths`.
  const int& nths = _nths;
  //! \brief Read only view of `_nph`.
  const int& nph = _nph;
  //! \brief Read only view of `_nphs`.
  const int& nphs = _nphs;
  //! \brief Read only view of `_nk`.
  const int& nk = _nk;
  //! \brief Read only view of `_nks`.
  const int& nks = _nks;
  //! \brief Read only view of `_nkks`.
  const int& nkks = _nkks;
  //! \brief Read only view of `_th`.
  const double& th = _th;
  //! \brief Read only view of `_thstp`.
  const double& thstp = _thstp;
  //! \brief Read only view of `_thlst`.
  const double& thlst = _thlst;
  //! \brief Read only view of `_ths`.
  const double& ths = _ths;
  //! \brief Read only view of `_thsstp`.
  const double& thsstp = _thsstp;
  //! \brief Read only view of `_thslst`.
  const double& thslst = _thslst;
  //! \brief Read only view of `_ph`.
  const double& ph = _ph;
  //! \brief Read only view of `_phstp`.
  const double& phstp = _phstp;
  //! \brief Read only view of `_phlst`.
  const double& phlst = _phlst;
  //! \brief Read only view of `_phs`.
  const double& phs = _phs;
  //! \brief Read only view of `_phsstp`.
  const double& phsstp = _phsstp;
  //! \brief Read only view of `_phslst`.
  const double& phslst = _phslst;
  //! \brief Read only view of `_thsca`.
  const double& thsca = _thsca;

  /*! \brief ScatteringAngles instance constructor.
   *
   * \param gconf: `GeometryConfiguration*` Pointer to a GeometryConfiguration object.
   */
  ScatteringAngles(GeometryConfiguration *gconf);
  
  /*! \brief ScatteringAngles copy constructor.
   *
   * \param rhs: `ScatteringAngles&` Reference to the ScatteringAngles object to be copied.
   */
  ScatteringAngles(const ScatteringAngles &rhs);

#ifdef MPI_VERSION
  /*! \brief ScatteringAngles copy from MPI broadcast.
   *
   * \param mpidata: `mixMPI *` Pointer to the mpidata instance used to copy the data.
   */
  ScatteringAngles(const mixMPI *mpidata);

    /*! \brief Broadcast over MPI the ScatteringAngles instance from MPI process 0 to all others.
   *
   * When using MPI, the initial ScatteringAngles instance created by MPI process 0
   * needs to be replicated on all other processes. This function sends it using
   * MPI broadcast calls. The MPI broadcast calls in this function must match those
   * in the constructor using the mixMPI pointer.
   *
   * \param mpidata: `mixMPI *` Pointer to the mpi structure used to do the MPI broadcast.
   */
  void mpibcast(const mixMPI *mpidata);
#endif

};


#endif
