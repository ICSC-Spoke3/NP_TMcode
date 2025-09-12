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

/*! \file outputs.h
 *
 * \brief Definition of the output format system.
 */
#ifndef INCLUDE_OUTPUTS_H_
#define INCLUDE_OUTPUTS_H_

// >>> OUTPUT FOR CLUSTER <<<
/*! \brief Class to collect output information for scattering from clusters.
 *
 * The results of the calculation can be saved in different formats.
 * It is therefore convenient to have a proper memory structure that
 * allows for storing the results and flushing them in any of the
 * permitted formats with just one operation. The purpose of the
 * `ClusterOutputInfo` class is to provide a wrapper for the output
 * of the cluster scattering solver.
 */
class ClusterOutputInfo {
protected:
  //! \brief Flag for skipping mpisend() and mpireceive()
  int _skip_flag;
  //! \brief Number of incident azimuth calculations.
  int _num_theta;
  //! \brief Number of scattered azimuth calculations.
  int _num_thetas;
  //! \brief Number of incident elevation calculations.
  int _num_phi;
  //! \brief Number of scattered elevation calculations.
  int _num_phis;
  //! \brief ID of the first computed wavelength
  int _first_xi;
  
  /*! \brief Write the output to a HDF5 file.
   *
   * \param file_name: `const string &` Path to the output to be written.
   * \return result: `int` Exit code (0 if successful).
   */
  int write_hdf5(const std::string &file_name);
  
  /*! \brief Write the output to a legacy text file.
   *
   * This function takes care of writing the output using the legacy
   * formatted ASCII structure. If the output file does not exist, it
   * is created. If it exists, the new content overwritten.
   *
   * \param output: `const string &` Path to the output to be written.
   * \return result: `int` Exit code (0 if successful).
   */
  int write_legacy(const std::string &output);
  
public:
  //! \brief Read-only view on skip_flag
  const int &skip_flag = _skip_flag;
  //! \brief Read-only view on the ID of the first scale
  const int &first_xi = _first_xi;
  //! \brief Number of spheres in the aggregate.
  int nsph;
  //! \brief Maximum internal field expansion order.
  int li;
  //! \brief Maximum external field expansion order.
  int le;
  //! \brief Maximum field expansion order.
  int lm;
  //! \brief Maximum coefficient matrix dimension.
  np_int mxndm;
  //! \brief Incident polarization flag.
  int inpol;
  //! \brief Number of points for transition layer integration.
  int npnt;
  //! \brief Number of points for non-transition layer integration.
  int npntts;
  //! \brief Flag for intensity.
  int iavm;
  //! \brief Flag for reference to meridional plane.
  int isam;
  //! \brief Flag for dielectric function definition.
  int idfc;
  //! \brief Vector of spherical components X Cartesian coordinates.
  double *vec_x_coords;
  //! \brief Vector of spherical components Y Cartesian coordinates.
  double *vec_y_coords;
  //! \brief Vector of spherical components Z Cartesian coordinates.
  double *vec_z_coords;
  //! \brief First incident radiation azimuth angle.
  double th;
  //! \brief Incident radiation azimuth angle step.
  double thstp;
  //! \brief Last incident radiation azimuth angle.
  double thlst;
  //! \brief First scattered radiation azimuth angle.
  double ths;
  //! \brief Scattered radiation azimuth angle step.
  double thsstp;
  //! \brief Last scattered radiation azimuth angle.
  double thslst;
  //! \brief First incident radiation elevation angle.
  double ph;
  //! \brief Incident radiation elevation angle step.
  double phstp;
  //! \brief Last incident radiation elevation angle.
  double phlst;
  //! \brief First scattered radiation elevation angle.
  double phs;
  //! \brief Scattered radiation elevation angle step.
  double phsstp;
  //! \brief Last scattered radiation elevation angle.
  double phslst;
  //! \brief Number of directions to be explicitly solved.
  int ndirs;
  //! \brief Refractive index of external medium.
  double exri;
  //! \brief Number of scales (wavelengths)
  int nxi;
  //! \brief Number of scales handled by the current process.
  int xi_block_size;
  //! \brief Index of the wavelength for T-matrix output.
  int jwtm;
  //! \brief Vector of scale (wavelength) indices.
  int *vec_jxi;
  //! \brief Vector of error severities (0 - success, 1 - HJV, 2 - DME).
  short *vec_ier;
  //! \brief Vector of vacuum wave numbers.
  double *vec_vk;
  //! \brief Vector of computed scales.
  double *vec_xi;
  //! \brief Number of sphere configurations.
  int configurations;
  //! \brief Vector of sphere sizes (all configurations for every scale).
  double *vec_sphere_sizes;
  //! \brief Vector of sphere refractive indices  (all configurations for every scale).
  dcomplex *vec_sphere_ref_indices;
  //! \brief Vector of sphere scattering cross-sections (all configurations for every scale).
  double *vec_sphere_scs;
  //! \brief Vector of sphere absorption cross-sections (all configurations for every scale).
  double *vec_sphere_abs;
  //! \brief Vector of sphere extinction cross-sections (all configurations for every scale).
  double *vec_sphere_exs;
  //! \brief Vector of sphere albedos (all configurations for every scale).
  double *vec_sphere_albs;
  //! \brief Vector of sphere scattering cross-section to geometric section ratios (all configurations for every scale).
  double *vec_sphere_sqscs;
  //! \brief Vector of sphere absorption cross-sections to geometric section ratios (all configurations for every scale).
  double *vec_sphere_sqabs;
  //! \brief Vector of sphere extinction cross-sections to geometric section ratios (all configurations for every scale).
  double *vec_sphere_sqexs;
  //! \brief Vector of sphere forward scattering amplitudes (all configurations for every scale).
  dcomplex *vec_fsas;
  //! \brief Vector of QSCHU = 4 pi IMAG(FSAS) / TOTAL_GEOM_SECTION (all configurations for every scale).
  double *vec_qschus;
  //! \brief Vector of PSCHU = 4 pi REAL(FSAS) / TOTAL_GEOM_SECTION (all configurations for every scale).
  double *vec_pschus;
  //! \brief Vector of S0MAG = ABS(FSAS) / (4 pi k^3) (all configurations for every scale).
  double *vec_s0mags;
  //! \brief Vector of sphere asymmetry parameters (all configurations for every scale).
  double *vec_cosavs;
  //! \brief Vector of sphere radiation pressure forces (all configurations for every scale).
  double *vec_raprs;
  //! \brief Vector of extinction contributions to radiation torques along k for parallel linear polarization (all configurations for every scale).
  double *vec_tqek1;
  //! \brief Vector of scattering contributions to radiation torques along k for parallel linear polarization (all configurations for every scale).
  double *vec_tqsk1;
  //! \brief Vector of extinction contributions to radiation torques along k for perpendicular linear polarization (all configurations for every scale).
  double *vec_tqek2;
  //! \brief Vector of scattering contributions to radiation torques along k for perpendicular linear polarization (all configurations for every scale).
  double *vec_tqsk2;
  //! \brief Vector of total forward scattering amplitudes (one for each scale).
  dcomplex *vec_fsat;
  //! \brief Vector of total QSCHU (one for each scale).
  double *vec_qschut;
  //! \brief Vector of total PSCHU (one for each scale).
  double *vec_pschut;
  //! \brief Vector of total S0MAG (one for each scale).
  double *vec_s0magt;
  //! \brief Total geometric section.
  double tgs;
  //! \brief Vector of cluster scattering cross-sections (parallel polarization).
  double *vec_scc1;
  //! \brief Vector of cluster scattering cross-sections (perpendicular polarization).
  double *vec_scc2;
  //! \brief Vector of cluster absorption cross-sections (parallel polarization).
  double *vec_abc1;
  //! \brief Vector of cluster absorption cross-sections (perpendicular polarization).
  double *vec_abc2;
  //! \brief Vector of cluster extinction cross-sections (parallel polarization).
  double *vec_exc1;
  //! \brief Vector of cluster extinction cross-sections (perpendicular polarization).
  double *vec_exc2;
  //! \brief Vector of cluster albedos (parallel polarization).
  double *vec_albedc1;
  //! \brief Vector of cluster albedos (perpendicular polarization).
  double *vec_albedc2;
  //! \brief Vector of cluster scattering to geometric cross-section ratios (parallel polarization).
  double *vec_qscamc1;
  //! \brief Vector of cluster scattering to geometric cross-section ratios (perpendicular polarization).
  double *vec_qscamc2;
  //! \brief Vector of cluster absorption to geometric cross-section ratios (parallel polarization).
  double *vec_qabsmc1;
  //! \brief Vector of cluster absorption to geometric cross-section ratios (perpendicular polarization).
  double *vec_qabsmc2;
  //! \brief Vector of cluster extinction to geometric cross-section ratios (parallel polarization).
  double *vec_qextmc1;
  //! \brief Vector of cluster extinction to geometric cross-section ratios (perpendicular polarization).
  double *vec_qextmc2;
  //! \brief Vector of cluster-to-sum-of-spheres scattering cross-section ratios (parallel polarization).
  double *vec_sccrt1;
  //! \brief Vector of cluster-to-sum-of-spheres scattering cross-section ratios (parallel polarization).
  double *vec_sccrt2;
  //! \brief Vector of cluster-to-sum-of-spheres absorption cross-section ratios (parallel polarization).
  double *vec_abcrt1;
  //! \brief Vector of cluster-to-sum-of-spheres absorption cross-section ratios (perpendicular polarization).
  double *vec_abcrt2;
  //! \brief Vector of cluster-to-sum-of-spheres extinction cross-section ratios (parallel polarization).
  double *vec_excrt1;
  //! \brief Vector of cluster-to-sum-of-spheres extinction cross-section ratios (perpendicular polarization).
  double *vec_excrt2;
  //! \brief Vector of forward scattering amplitudes for polarization parallel to incidence (one per scale).
  dcomplex *vec_fsac11;
  //! \brief Vector of forward scattering amplitudes for polarization perpendicular to incidence (one per scale).
  dcomplex *vec_fsac21;
  //! \brief Vector of forward scattering amplitudes for polarization parallel to incidence (one per scale).
  dcomplex *vec_fsac22;
  //! \brief Vector of forward scattering amplitudes for polarization perpendicular to incidence (one per scale).
  dcomplex *vec_fsac12;
  //! \brief Vector of cluster QSCHU (parallel polarization).
  double *vec_qschuc1;
  //! \brief Vector of cluster QSCHU (perpendicular polarization).
  double *vec_qschuc2;
  //! \brief Vector of cluster PSCHU (parallel polarization).
  double *vec_pschuc1;
  //! \brief Vector of cluster PSCHU (perpendicular polarization).
  double *vec_pschuc2;
  //! \brief Vector of cluster S0MAG (parallel polarization).
  double *vec_s0magc1;
  //! \brief Vector of cluster S0MAG (perpendicular polarization).
  double *vec_s0magc2;
  //! \brief Vector of cluster asymmetry parameters (parallel polarization).
  double *vec_cosavc1;
  //! \brief Vector of cluster asymmetry parameters (perpendicular polarization).
  double *vec_cosavc2;
  //! \brief Vector of cluster radiation pressure forces (parallel polarization).
  double *vec_raprc1;
  //! \brief Vector of cluster radiation pressure forces (perpendicular polarization).
  double *vec_raprc2;
  //! \brief Vector of optical forces along incidence direction [N] (parallel polarization).
  double *vec_fkc1;
  //! \brief Vector of optical forces along incidence direction [N] (perpendicular polarization).
  double *vec_fkc2;
  //! \brief Vector of incidence azimuth directions (one per incidence azimuth).
  double *vec_dir_tidg;
  //! \brief Vector of incidence elevation directions (one per incidence elevation).
  double *vec_dir_pidg;
  //! \brief Vector of scattering azimuth directions (one per scattering azimuth).
  double *vec_dir_tsdg;
  //! \brief Vector of scattering elevation directions (one per scattering elevation).
  double *vec_dir_psdg;
  //! \brief Vector of scattering angles (one per direction).
  double *vec_dir_scand;
  //! \brief Control parameter for incidence plane referred to meridional plane (one per direction).
  double *vec_dir_cfmp;
  //! \brief Control parameter for scattering plane referred to meridional plane (one per direction).
  double *vec_dir_sfmp;
  //! \brief Control parameter for incidence plane referred to scattering plane (one per direction).
  double *vec_dir_cfsp;
  //! \brief Control parameter for scattering plane referred to scattering plane (one per direction).
  double *vec_dir_sfsp;
  //! \brief Components of the unitary vector perpendicular to incidence plane (three per direction).
  double *vec_dir_un;
  //! \brief Components of the unitary vector perpendicular to scattering plane (three per direction).
  double *vec_dir_uns;
  //! \brief Vector of sphere differential scattering amplitude with polarization parallel to parallel incidence field.
  dcomplex *vec_dir_sas11;
  //! \brief Vector of sphere differential scattering amplitude with polarization perpendicular to the parallel incidence field.
  dcomplex *vec_dir_sas21;
  //! \brief Vector of sphere differential scattering amplitude with polarization perpendicular to perpendicular incidence field.
  dcomplex *vec_dir_sas12;
  //! \brief Vector of sphere differential scattering amplitude with polarization parallel the perpendicular incidence field.
  dcomplex *vec_dir_sas22;
  //! \brief Vector of sphere Mueller transormation matrices referred to meridional plane.
  double *vec_dir_muls;
  //! \brief Vector of sphere Mueller transormation matrices referred to scattering plane.
  double *vec_dir_mulslr;
  //! \brief Vector of sphere total differential scattering amplitude with polarization parallel to parallel incidence field.
  dcomplex *vec_dir_sat11;
  //! \brief Vector of sphere total differential scattering amplitude with polarization perpendicular to the parallel incidence field.
  dcomplex *vec_dir_sat21;
  //! \brief Vector of sphere total differential scattering amplitude with polarization perpendicular to perpendicular incidence field.
  dcomplex *vec_dir_sat12;
  //! \brief Vector of sphere total differential scattering amplitude with polarization parallel the perpendicular incidence field.
  dcomplex *vec_dir_sat22;
  //! \brief Vector of cluster differential scattering cross-sections (parallel polarization).
  double *vec_dir_scc1;
  //! \brief Vector of cluster differential scattering cross-sections (perpendicular polarization).
  double *vec_dir_scc2;
  //! \brief Vector of cluster differential absorption cross-sections (parallel polarization).
  double *vec_dir_abc1;
  //! \brief Vector of cluster differential absorption cross-sections (perpendicular polarization).
  double *vec_dir_abc2;
  //! \brief Vector of cluster differential extinction cross-sections (parallel polarization).
  double *vec_dir_exc1;
  //! \brief Vector of cluster differential extinction cross-sections (perpendicular polarization).
  double *vec_dir_exc2;
  //! \brief Vector of cluster differential albedos (parallel polarization).
  double *vec_dir_albedc1;
  //! \brief Vector of cluster differential albedos (perpendicular polarization).
  double *vec_dir_albedc2;
  //! \brief Vector of differential scattering to geometric cross-section ratios (parallel polarization).
  double *vec_dir_qscc1;
  //! \brief Vector of differential scattering to geometric cross-section ratios (perpendicular polarization).
  double *vec_dir_qscc2;
  //! \brief Vector of differential absorption to geometric cross-section ratios (parallel polarization).
  double *vec_dir_qabc1;
  //! \brief Vector of differential absorption to geometric cross-section ratios (perpendicular polarization).
  double *vec_dir_qabc2;
  //! \brief Vector of differential extinction to geometric cross-section ratios (parallel polarization).
  double *vec_dir_qexc1;
  //! \brief Vector of differential extinction to geometric cross-section ratios (perpendicular polarization).
  double *vec_dir_qexc2;
  //! \brief Vector of differential cluster-to-total scattering cross-section ratios (parallel polarization).
  double *vec_dir_sccrt1;
  //! \brief Vector of differential cluster-to-total scattering cross-section ratios (perpendicular polarization).
  double *vec_dir_sccrt2;
  //! \brief Vector of differential cluster-to-total absorption cross-section ratios (parallel polarization).
  double *vec_dir_abcrt1;
  //! \brief Vector of differential cluster-to-total absorption cross-section ratios (perpendicular polarization).
  double *vec_dir_abcrt2;
  //! \brief Vector of differential cluster-to-total extinction cross-section ratios (parallel polarization).
  double *vec_dir_excrt1;
  //! \brief Vector of differential cluster-to-total extinction cross-section ratios (perpendicular polarization).
  double *vec_dir_excrt2;
  //! \brief Vector of differential cluster forward scattering amplitude with polarization parallel to parallel incidence field (one per direction and scale).
  dcomplex *vec_dir_fsac11;
  //! \brief Vector of differential cluster forward scattering amplitude with polarization perpendicular to the parallel incidence field (one per direction and scale).
  dcomplex *vec_dir_fsac21;
  //! \brief Vector of differential cluster forward scattering amplitude with polarization perpendicular to perpendicular incidence field (one per direction and scale).
  dcomplex *vec_dir_fsac12;
  //! \brief Vector of differential cluster forward scattering amplitude with polarization parallel the perpendicular incidence field (one per direction and scale).
  dcomplex *vec_dir_fsac22;
  //! \brief Vector of differential cluster scattering amplitude with polarization parallel to parallel incidence field (one per direction and scale).
  dcomplex *vec_dir_sac11;
  //! \brief Vector of differential cluster scattering amplitude with polarization perpendicular to the parallel incidence field (one per direction and scale).
  dcomplex *vec_dir_sac21;
  //! \brief Vector of differential cluster scattering amplitude with polarization perpendicular to perpendicular incidence field (one per direction and scale).
  dcomplex *vec_dir_sac12;
  //! \brief Vector of differential cluster scattering amplitude with polarization parallel the perpendicular incidence field (one per direction and scale).
  dcomplex *vec_dir_sac22;
  //! \brief Vector of differential cluster QSCHU (parallel polarization).
  double *vec_dir_qschuc1;
  //! \brief Vector of differential cluster QSCHU (perpendicular polarization).
  double *vec_dir_qschuc2;
  //! \brief Vector of differential cluster PSCHU (parallel polarization).
  double *vec_dir_pschuc1;
  //! \brief Vector of differential cluster PSCHU (perpendicular polarization).
  double *vec_dir_pschuc2;
  //! \brief Vector of cluster differential S0MAG (parallel polarization).
  double *vec_dir_s0magc1;
  //! \brief Vector of cluster differential S0MAG (perpendicular polarization).
  double *vec_dir_s0magc2;
  //! \brief Vector of differential cluster asymmetry parameters (parallel polarization).
  double *vec_dir_cosavc1;
  //! \brief Vector of differential cluster asymmetry parameters (perpendicular polarization).
  double *vec_dir_cosavc2;
  //! \brief Vector of differential cluster radiation pressure forces (1).
  double *vec_dir_raprc1;
  //! \brief Vector of differential cluster radiation pressure forces (1).
  double *vec_dir_raprc2;
  //! \brief Vector of differential radiation pressure force components along the polarization direction (parallel polarization).
  double *vec_dir_flc1;
  //! \brief Vector of differential radiation pressure force components along the polarization direction (perpendicular polarization).
  double *vec_dir_flc2;
  //! \brief Vector of differential radiation pressure force components perpendicular to the polarization direction (parallel polarization).
  double *vec_dir_frc1;
  //! \brief Vector of differential radiation pressure force components perpendicular to the polarization direction (perpendicular polarization).
  double *vec_dir_frc2;
  //! \brief Vector of differential radiation pressure force components along the incidence direction (parallel polarization).
  double *vec_dir_fkc1;
  //! \brief Vector of differential radiation pressure force components along the incidence direction (perpendicular polarization).
  double *vec_dir_fkc2;
  //! \brief Vector of differential radiation pressure force components along the X axis (parallel polarization).
  double *vec_dir_fxc1;
  //! \brief Vector of differential radiation pressure force components along the X axis (perpendicular polarization).
  double *vec_dir_fxc2;
  //! \brief Vector of differential radiation pressure force components along the Y axis (parallel polarization).
  double *vec_dir_fyc1;
  //! \brief Vector of differential radiation pressure force components along the Y axis (perpendicular polarization).
  double *vec_dir_fyc2;
  //! \brief Vector of differential radiation pressure force components along the Z axis (parallel polarization).
  double *vec_dir_fzc1;
  //! \brief Vector of differential radiation pressure force components along the Z axis (perpendicular polarization).
  double *vec_dir_fzc2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the polarization direction (parallel polarization).
  double *vec_dir_tqelc1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the polarization direction (perpendicular polarization).
  double *vec_dir_tqelc2;
  //! \brief Vector of differential extinction contribution to radiation torque components perpendicular to the polarization direction (parallel polarization).
  double *vec_dir_tqerc1;
  //! \brief Vector of differential extinction contribution to radiation torque components perpendicular to the polarization direction (perpendicular polarization).
  double *vec_dir_tqerc2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the incidence direction (parallel polarization).
  double *vec_dir_tqekc1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the incidence direction (perpendicular polarization).
  double *vec_dir_tqekc2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the X axis (parallel polarization).
  double *vec_dir_tqexc1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the X axis (perpendicular polarization).
  double *vec_dir_tqexc2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Y axis (parallel polarization).
  double *vec_dir_tqeyc1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Y axis (perpendicular polarization).
  double *vec_dir_tqeyc2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Z axis (parallel polarization).
  double *vec_dir_tqezc1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Z axis (perpendicular polarization).
  double *vec_dir_tqezc2;
  //! \brief Vector of differential scattering contribution to radiation torque components along the polarization direction (parallel polarization).
  double *vec_dir_tqslc1;
  //! \brief Vector of differential scattering contribution to radiation torque components along the polarization direction (perpendicular polarization).
  double *vec_dir_tqslc2;
  //! \brief Vector of differential scattering contribution to radiation torque components perpendicular to the polarization direction (parallel polarization).
  double *vec_dir_tqsrc1;
  //! \brief Vector of differential scattering contribution to radiation torque components perpendicular to the polarization direction (perpendicular polarization).
  double *vec_dir_tqsrc2;
  //! \brief Vector of differential scattering contribution to radiation torque components along the incidence direction (parallel polarization).
  double *vec_dir_tqskc1;
  //! \brief Vector of differential scattering contribution to radiation torque components along the incidence direction (perpendicular polarization).
  double *vec_dir_tqskc2;
  //! \brief Vector of differential scattering contribution to radiation torque components along X axis (parallel polarization).
  double *vec_dir_tqsxc1;
  //! \brief Vector of differential scattering contribution to radiation torque components along X axis (perpendicular polarization).
  double *vec_dir_tqsxc2;
  //! \brief Vector of differential scattering contribution to radiation torque components along Y axis (parallel polarization).
  double *vec_dir_tqsyc1;
  //! \brief Vector of differential scattering contribution to radiation torque components along Y axis (perpendicular polarization).
  double *vec_dir_tqsyc2;
  //! \brief Vector of differential scattering contribution to radiation torque components along Z axis (parallel polarization).
  double *vec_dir_tqszc1;
  //! \brief Vector of differential scattering contribution to radiation torque components along Z axis (perpendicular polarization).
  double *vec_dir_tqszc2;
  //! \brief Vector of cluster Mueller transormation matrices referred to meridional plane (16 per direction per scale).
  double *vec_dir_mulc;
  //! \brief Vector of cluster Mueller transormation matrices referred to scattering plane (16 per direction per scale).
  double *vec_dir_mulclr;
  
  /*! \brief `ClusterOutputInfo` default instance constructor.
   *
   * \param sc: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` instance.
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \param mpidata: `const mixMPI*` Pointer to a mixMPI instance.
   * \param first_xi: `int` Index of the first scale in output (optional, default is 1).
   * \param xi_length: `int` Number of scales tobe included in output (optional, default is 0, meaning all).
   */   
  ClusterOutputInfo(
    ScattererConfiguration *sc, GeometryConfiguration *gc,
    const mixMPI *mpidata, int first_xi = 1, int xi_length = 0
  );

  /*! \brief `ClusterOutputInfo` constructor from HDF5 input.
   *
   * \param hdf5_name: `const string &` Path to the HDF5 file to be read.
   */   
  ClusterOutputInfo(const std::string &hdf5_name);

  /*! \brief `ClusterOutputInfo` constructor for the dummy NULL case
   *
   * \param skip_flag: `const int` must be passed as the `1` constant.
   */   
  ClusterOutputInfo(const int skip_flag);

  /*! \brief `ClusterOutputInfo` instance destroyer.
   */
  ~ClusterOutputInfo();

  /*! \brief Estimate the size of the structure that would be built for given input.
   *
   * \param sc: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` instance.
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \param first_xi: `int` Index of the first scale in output (optional, default is 1).
   * \param xi_length: `int` Number of scales tobe included in output (optional, default is all).
   * \return size: `long` Estimated instance size in bytes.
   */
  static long compute_size(
    ScattererConfiguration *sc, GeometryConfiguration *gc,
    int first_xi = 1, int xi_length = 0
  );
  
  /*! \brief Get the size of a `ClusterOutputInfo` instance in bytes.
   *
   * \return size: `long` Estimated instance size in bytes.
   */
  long compute_size();

  /*! \brief Insert in the current output data the data of another block.
   *
   * \param rhs: `const ClusterOutputInfo &` Reference to the source data block.
   * \return result: `int` Exit code (0 if successful).
   */
  int insert(const ClusterOutputInfo &rhs);

  /*! \brief Write the output to a file.
   *
   * \param output: `const string &` Path to the output to be written.
   * \param format: `const string &` Output format (one of LEGACY or HDF5).
   * \return result: `int` Exit code (0 if successful).
   */
  int write(const std::string &output, const std::string &format);

#ifdef MPI_VERSION
  /*! \brief Receive output data from worker processes.
   *
   * This function is invoked by the MPI rank-0 process to fetch the
   * output data produced by higher rank processes. When calling this
   * function, process 0 halts until a valid data chunk is transmitted
   * by the queried process.
   *
   * \param mpidata: `const mixMPI*` Pointer to a `mixMPI` instance.
   * \param pid: `int` Rank of the process that is transmitting data.
   * \return result: `int` An exit code (0 for success).
   */
  int mpireceive(const mixMPI* mpidata, int pid);

  /*! \brief Send output data to process 0.
   *
   * This function is invoked by non-zero ranked MPI processes when
   * they are ready to send back the output data. When a process meets
   * this function call, it halts until MPI process 0 asks for the
   * data transmission.
   *
   * \param mpidata: `const mixMPI*` Pointer to a `mixMPI` instance.
   * \param pid: `int` Rank of the process that is transmitting data.
   * \return result: `int` An exit code (0 for success).
   */
  int mpisend(const mixMPI *mpidata);
#endif // MPI_VERSION
};
// >>> END OF OUTPUT FOR CLUSTER <<<

// >>> OUTPUT FOR INCLUSION <<<
/*! \brief Class to collect output information for scattering from particle with inclusions.
 *
 * The results of the calculation can be saved in different formats.
 * It is therefore convenient to have a proper memory structure that
 * allows for storing the results and flushing them in any of the
 * permitted formats with just one operation. The purpose of the
 * `InclusionOutputInfo` class is to provide a wrapper for the output
 * of the particle with inclusions scattering solver.
 */
class InclusionOutputInfo {
protected:
  //! \brief Flag for skipping mpisend() and mpireceive()
  int _skip_flag;
  //! \brief Number of incident azimuth calculations.
  int _num_theta;
  //! \brief Number of scattered azimuth calculations.
  int _num_thetas;
  //! \brief Number of incident elevation calculations.
  int _num_phi;
  //! \brief Number of scattered elevation calculations.
  int _num_phis;
  //! \brief ID of the first computed wavelength
  int _first_xi;
  
  /*! \brief Write the output to a HDF5 file.
   *
   * \param file_name: `const string &` Path to the output to be written.
   * \return result: `int` Exit code (0 if successful).
   */
  int write_hdf5(const std::string &file_name);
  
  /*! \brief Write the output to a legacy text file.
   *
   * This function takes care of writing the output using the legacy
   * formatted ASCII structure. If the output file does not exist, it
   * is created. If it exists, the new content is overwritten.
   *
   * \param output: `const string &` Path to the output to be written.
   * \return result: `int` Exit code (0 if successful).
   */
  int write_legacy(const std::string &output);
  
public:
  //! \brief Read-only view on skip_flag
  const int &skip_flag = _skip_flag;
  //! \brief Read-only view on the ID of the first scale
  const int &first_xi = _first_xi;
  //! \brief Number of spheres in the aggregate.
  int nsph;
  //! \brief Maximum internal field expansion order.
  int li;
  //! \brief Maximum external field expansion order.
  int le;
  //! \brief Maximum field expansion order.
  int lm;
  //! \brief Maximum coefficient matrix dimension.
  np_int mxndm;
  //! \brief Incident polarization flag.
  int inpol;
  //! \brief Number of points for transition layer integration.
  int npnt;
  //! \brief Number of points for non-transition layer integration.
  int npntts;
  //! \brief Flag for intensity.
  int iavm;
  //! \brief Flag for reference to meridional plane.
  int isam;
  //! \brief Flag for dielectric function definition.
  int idfc;
  //! \brief Vector of spherical components X Cartesian coordinates.
  double *vec_x_coords;
  //! \brief Vector of spherical components Y Cartesian coordinates.
  double *vec_y_coords;
  //! \brief Vector of spherical components Z Cartesian coordinates.
  double *vec_z_coords;
  //! \brief First incident radiation azimuth angle.
  double th;
  //! \brief Incident radiation azimuth angle step.
  double thstp;
  //! \brief Last incident radiation azimuth angle.
  double thlst;
  //! \brief First scattered radiation azimuth angle.
  double ths;
  //! \brief Scattered radiation azimuth angle step.
  double thsstp;
  //! \brief Last scattered radiation azimuth angle.
  double thslst;
  //! \brief First incident radiation elevation angle.
  double ph;
  //! \brief Incident radiation elevation angle step.
  double phstp;
  //! \brief Last incident radiation elevation angle.
  double phlst;
  //! \brief First scattered radiation elevation angle.
  double phs;
  //! \brief Scattered radiation elevation angle step.
  double phsstp;
  //! \brief Last scattered radiation elevation angle.
  double phslst;
  //! \brief Number of directions to be explicitly solved.
  int ndirs;
  //! \brief Refractive index of external medium.
  double exri;
  //! \brief Number of scales (wavelengths)
  int nxi;
  //! \brief Number of scales handled by the current process.
  int xi_block_size;
  //! \brief Index of the wavelength for T-matrix output.
  int jwtm;
  //! \brief Vector of scale (wavelength) indices.
  int *vec_jxi;
  //! \brief Vector of error severities (0 - success, 1 - INDME, 2 - OSPV).
  short *vec_ier;
  //! \brief Vector of vacuum wave numbers.
  double *vec_vk;
  //! \brief Vector of computed scales.
  double *vec_xi;
  //! \brief Number of sphere configurations.
  int configurations;
  //! \brief Vector of sphere sizes (all configurations for every scale).
  double *vec_sphere_sizes;
  //! \brief Vector of sphere refractive indices  (all configurations for every scale).
  dcomplex *vec_sphere_ref_indices;
  //! \brief Vector of particle scattering cross-sections (parallel polarization).
  double *vec_scs1;
  //! \brief Vector of particle scattering cross-sections (perpendicular polarization).
  double *vec_scs2;
  //! \brief Vector of particle absorption cross-sections (parallel polarization).
  double *vec_abs1;
  //! \brief Vector of particle absorption cross-sections (perpendicular polarization).
  double *vec_abs2;
  //! \brief Vector of particle extinction cross-sections (parallel polarization).
  double *vec_exs1;
  //! \brief Vector of particle extinction cross-sections (perpendicular polarization).
  double *vec_exs2;
  //! \brief Vector of particle albedos (parallel polarization).
  double *vec_albeds1;
  //! \brief Vector of particle albedos (perpendicular polarization).
  double *vec_albeds2;
  //! \brief Vector of particle scattering-to-geometric cross-sections (parallel polarization).
  double *vec_scsrt1;
  //! \brief Vector of particle scattering-to-geometric cross-sections (perpendicular polarization).
  double *vec_scsrt2;
  //! \brief Vector of particle absorption-to-geometric cross-sections (parallel polarization).
  double *vec_absrt1;
  //! \brief Vector of particle absorption-to-geometric cross-sections (perpendicular polarization).
  double *vec_absrt2;
  //! \brief Vector of particle extinction-to-geometric cross-sections (parallel polarization).
  double *vec_exsrt1;
  //! \brief Vector of particle extinction-to-geometric cross-sections (perpendicular polarization).
  double *vec_exsrt2;
  //! \brief Vector of particle QSCHU (parallel polarization).
  double *vec_qschu1;
  //! \brief Vector of particle QSCHU (perpendicular polarization).
  double *vec_qschu2;
  //! \brief Vector of particle PSCHU (parallel polarization).
  double *vec_pschu1;
  //! \brief Vector of particle PSCHU (perpendicular polarization).
  double *vec_pschu2;
  //! \brief Vector of particle S0MAG (parallel polarization).
  double *vec_s0mag1;
  //! \brief Vector of particle S0MAG (perpendicular polarization).
  double *vec_s0mag2;
  //! \brief Vector of particle average asymmetry parameter (parallel polarization).
  double *vec_cosav1;
  //! \brief Vector of particle average asymmetry parameter (perpendicular polarization).
  double *vec_cosav2;
  //! \brief Vector of particle average radiation pressure force (N - parallel polarization).
  double *vec_raprs1;
  //! \brief Vector of particle average radiation pressure force (N - perpendicular polarization).
  double *vec_raprs2;
  //! \brief Vector of particle average radiation force along incidence direction (N - parallel polarization).
  double *vec_fk1;
  //! \brief Vector of particle average radiation force along incidence direction (N - perpendicular polarization).
  double *vec_fk2;
  //! \brief Vector of forward scattering amplitudes for polarization parallel to incidence (one per scale).
  dcomplex *vec_fsas11;
  //! \brief Vector of forward scattering amplitudes for polarization perpendicular to incidence (one per scale).
  dcomplex *vec_fsas21;
  //! \brief Vector of forward scattering amplitudes for polarization parallel to incidence (one per scale).
  dcomplex *vec_fsas22;
  //! \brief Vector of forward scattering amplitudes for polarization perpendicular to incidence (one per scale).
  dcomplex *vec_fsas12;
  //! \brief Vector of incidence azimuth directions (one per incidence azimuth).
  double *vec_dir_tidg;
  //! \brief Vector of incidence elevation directions (one per incidence elevation).
  double *vec_dir_pidg;
  //! \brief Vector of scattering azimuth directions (one per scattering azimuth).
  double *vec_dir_tsdg;
  //! \brief Vector of scattering elevation directions (one per scattering elevation).
  double *vec_dir_psdg;
  //! \brief Vector of scattering angles (one per direction).
  double *vec_dir_scand;
  //! \brief Control parameter for incidence plane referred to meridional plane (one per direction).
  double *vec_dir_cfmp;
  //! \brief Control parameter for scattering plane referred to meridional plane (one per direction).
  double *vec_dir_sfmp;
  //! \brief Control parameter for incidence plane referred to scattering plane (one per direction).
  double *vec_dir_cfsp;
  //! \brief Control parameter for scattering plane referred to scattering plane (one per direction).
  double *vec_dir_sfsp;
  //! \brief Components of the unitary vector perpendicular to incidence plane (three per direction).
  double *vec_dir_un;
  //! \brief Components of the unitary vector perpendicular to scattering plane (three per direction).
  double *vec_dir_uns;
  //! \brief Vector of particle differential scattering cross-sections (parallel polarization).
  double *vec_dir_scs1;
  //! \brief Vector of particle differential scattering cross-sections (perpendicular polarization).
  double *vec_dir_scs2;
  //! \brief Vector of particle differential absorption cross-sections (parallel polarization).
  double *vec_dir_abs1;
  //! \brief Vector of particle differential absorption cross-sections (perpendicular polarization).
  double *vec_dir_abs2;
  //! \brief Vector of particle differential extinction cross-sections (parallel polarization).
  double *vec_dir_exs1;
  //! \brief Vector of particle differential extinction cross-sections (perpendicular polarization).
  double *vec_dir_exs2;
  //! \brief Vector of particle differential albedos (parallel polarization).
  double *vec_dir_albeds1;
  //! \brief Vector of particle differential albedos (perpendicular polarization).
  double *vec_dir_albeds2;
  //! \brief Vector of particle differential scattering-to-geometric cross-sections (parallel polarization).
  double *vec_dir_scsrt1;
  //! \brief Vector of particle differential scattering-to-geometric cross-sections (perpendicular polarization).
  double *vec_dir_scsrt2;
  //! \brief Vector of particle differential absorption-to-geometric cross-sections (parallel polarization).
  double *vec_dir_absrt1;
  //! \brief Vector of particle differential absorption-to-geometric cross-sections (perpendicular polarization).
  double *vec_dir_absrt2;
  //! \brief Vector of particle differential extinction-to-geometric cross-sections (parallel polarization).
  double *vec_dir_exsrt1;
  //! \brief Vector of particle differential extinction-to-geometric cross-sections (perpendicular polarization).
  double *vec_dir_exsrt2;
  //! \brief Vector of particle differential forward scattering amplitude with polarization parallel to parallel incidence field.
  dcomplex *vec_dir_fsas11;
  //! \brief Vector of particle differential forward scattering amplitude with polarization perpendicular to the parallel incidence field.
  dcomplex *vec_dir_fsas21;
  //! \brief Vector of particle differential forward scattering amplitude with polarization perpendicular to perpendicular incidence field.
  dcomplex *vec_dir_fsas12;
  //! \brief Vector of particle differential forward scattering amplitude with polarization parallel the perpendicular incidence field.
  dcomplex *vec_dir_fsas22;
  //! \brief Vector of particle differential scattering amplitude with polarization parallel to parallel incidence field.
  dcomplex *vec_dir_sas11;
  //! \brief Vector of particle differential scattering amplitude with polarization perpendicular to the parallel incidence field.
  dcomplex *vec_dir_sas21;
  //! \brief Vector of particle differential scattering amplitude with polarization perpendicular to perpendicular incidence field.
  dcomplex *vec_dir_sas12;
  //! \brief Vector of particle differential scattering amplitude with polarization parallel the perpendicular incidence field.
  dcomplex *vec_dir_sas22;
  //! \brief Vector of differential particle QSCHU (parallel polarization).
  double *vec_dir_qschu1;
  //! \brief Vector of differential particle QSCHU (perpendicular polarization).
  double *vec_dir_qschu2;
  //! \brief Vector of differential particle PSCHU (parallel polarization).
  double *vec_dir_pschu1;
  //! \brief Vector of differential particle PSCHU (perpendicular polarization).
  double *vec_dir_pschu2;
  //! \brief Vector of particle differential S0MAG (parallel polarization).
  double *vec_dir_s0mag1;
  //! \brief Vector of particle differential S0MAG (perpendicular polarization).
  double *vec_dir_s0mag2;
  //! \brief Vector of differential particle asymmetry parameters (parallel polarization).
  double *vec_dir_cosav1;
  //! \brief Vector of differential particle asymmetry parameters (perpendicular polarization).
  double *vec_dir_cosav2;
  //! \brief Vector of differential particle radiation pressure forces (1).
  double *vec_dir_rapr1;
  //! \brief Vector of differential particle radiation pressure forces (1).
  double *vec_dir_rapr2;
  //! \brief Vector of differential radiation pressure force components along the polarization direction (parallel polarization).
  double *vec_dir_fl1;
  //! \brief Vector of differential radiation pressure force components along the polarization direction (perpendicular polarization).
  double *vec_dir_fl2;
  //! \brief Vector of differential radiation pressure force components perpendicular to the polarization direction (parallel polarization).
  double *vec_dir_fr1;
  //! \brief Vector of differential radiation pressure force components perpendicular to the polarization direction (perpendicular polarization).
  double *vec_dir_fr2;
  //! \brief Vector of differential radiation pressure force components along the incidence direction (parallel polarization).
  double *vec_dir_fk1;
  //! \brief Vector of differential radiation pressure force components along the incidence direction (perpendicular polarization).
  double *vec_dir_fk2;
  //! \brief Vector of differential radiation pressure force components along the X axis (parallel polarization).
  double *vec_dir_fx1;
  //! \brief Vector of differential radiation pressure force components along the X axis (perpendicular polarization).
  double *vec_dir_fx2;
  //! \brief Vector of differential radiation pressure force components along the Y axis (parallel polarization).
  double *vec_dir_fy1;
  //! \brief Vector of differential radiation pressure force components along the Y axis (perpendicular polarization).
  double *vec_dir_fy2;
  //! \brief Vector of differential radiation pressure force components along the Z axis (parallel polarization).
  double *vec_dir_fz1;
  //! \brief Vector of differential radiation pressure force components along the Z axis (perpendicular polarization).
  double *vec_dir_fz2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the polarization direction (parallel polarization).
  double *vec_dir_tqel1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the polarization direction (perpendicular polarization).
  double *vec_dir_tqel2;
  //! \brief Vector of differential extinction contribution to radiation torque components perpendicular to the polarization direction (parallel polarization).
  double *vec_dir_tqer1;
  //! \brief Vector of differential extinction contribution to radiation torque components perpendicular to the polarization direction (perpendicular polarization).
  double *vec_dir_tqer2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the incidence direction (parallel polarization).
  double *vec_dir_tqek1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the incidence direction (perpendicular polarization).
  double *vec_dir_tqek2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the X axis (parallel polarization).
  double *vec_dir_tqex1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the X axis (perpendicular polarization).
  double *vec_dir_tqex2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Y axis (parallel polarization).
  double *vec_dir_tqey1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Y axis (perpendicular polarization).
  double *vec_dir_tqey2;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Z axis (parallel polarization).
  double *vec_dir_tqez1;
  //! \brief Vector of differential extinction contribution to radiation torque components along the Z axis (perpendicular polarization).
  double *vec_dir_tqez2;
  //! \brief Vector of differential scattering contribution to radiation torque components along the polarization direction (parallel polarization).
  double *vec_dir_tqsl1;
  //! \brief Vector of differential scattering contribution to radiation torque components along the polarization direction (perpendicular polarization).
  double *vec_dir_tqsl2;
  //! \brief Vector of differential scattering contribution to radiation torque components perpendicular to the polarization direction (parallel polarization).
  double *vec_dir_tqsr1;
  //! \brief Vector of differential scattering contribution to radiation torque components perpendicular to the polarization direction (perpendicular polarization).
  double *vec_dir_tqsr2;
  //! \brief Vector of differential scattering contribution to radiation torque components along the incidence direction (parallel polarization).
  double *vec_dir_tqsk1;
  //! \brief Vector of differential scattering contribution to radiation torque components along the incidence direction (perpendicular polarization).
  double *vec_dir_tqsk2;
  //! \brief Vector of differential scattering contribution to radiation torque components along X axis (parallel polarization).
  double *vec_dir_tqsx1;
  //! \brief Vector of differential scattering contribution to radiation torque components along X axis (perpendicular polarization).
  double *vec_dir_tqsx2;
  //! \brief Vector of differential scattering contribution to radiation torque components along Y axis (parallel polarization).
  double *vec_dir_tqsy1;
  //! \brief Vector of differential scattering contribution to radiation torque components along Y axis (perpendicular polarization).
  double *vec_dir_tqsy2;
  //! \brief Vector of differential scattering contribution to radiation torque components along Z axis (parallel polarization).
  double *vec_dir_tqsz1;
  //! \brief Vector of differential scattering contribution to radiation torque components along Z axis (perpendicular polarization).
  double *vec_dir_tqsz2;
  //! \brief Vector of cluster Mueller transormation matrices referred to meridional plane (16 per direction per scale).
  double *vec_dir_mull;
  //! \brief Vector of cluster Mueller transormation matrices referred to scattering plane (16 per direction per scale).
  double *vec_dir_mulllr;
  
  /*! \brief `InclusionOutputInfo` default instance constructor.
   *
   * \param sc: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` instance.
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \param mpidata: `const mixMPI*` Pointer to a mixMPI instance.
   * \param first_xi: `int` Index of the first scale in output (optional, default is 1).
   * \param xi_length: `int` Number of scales tobe included in output (optional, default is 0, meaning all).
   */   
  InclusionOutputInfo(
    ScattererConfiguration *sc, GeometryConfiguration *gc,
    const mixMPI *mpidata, int first_xi = 1, int xi_length = 0
  );

  /*! \brief `InclusionOutputInfo` constructor from HDF5 input.
   *
   * \param hdf5_name: `const string &` Path to the HDF5 file to be read.
   */   
  InclusionOutputInfo(const std::string &hdf5_name);

  /*! \brief `InclusionOutputInfo` constructor for the dummy NULL case
   *
   * \param skip_flag: `const int` must be passed as the `1` constant.
   */   
  InclusionOutputInfo(const int skip_flag);

  /*! \brief `InclusionOutputInfo` instance destroyer.
   */
  ~InclusionOutputInfo();

  /*! \brief Estimate the size of the structure that would be built for given input.
   *
   * \param sc: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` instance.
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \param first_xi: `int` Index of the first scale in output (optional, default is 1).
   * \param xi_length: `int` Number of scales tobe included in output (optional, default is all).
   * \return size: `long` Estimated instance size in bytes.
   */
  static long compute_size(
    ScattererConfiguration *sc, GeometryConfiguration *gc,
    int first_xi = 1, int xi_length = 0
  );
  
  /*! \brief Get the size of a `ClusterOutputInfo` instance in bytes.
   *
   * \return size: `long` Estimated instance size in bytes.
   */
  long compute_size();

  /*! \brief Insert in the current output data the data of another block.
   *
   * \param rhs: `const InclusionOutputInfo &` Reference to the source data block.
   * \return result: `int` Exit code (0 if successful).
   */
  int insert(const InclusionOutputInfo &rhs);

  /*! \brief Write the output to a file.
   *
   * \param output: `const string &` Path to the output to be written.
   * \param format: `const string &` Output format (one of LEGACY or HDF5).
   * \return result: `int` Exit code (0 if successful).
   */
  int write(const std::string &output, const std::string &format);

#ifdef MPI_VERSION
  /*! \brief Receive output data from worker processes.
   *
   * This function is invoked by the MPI rank-0 process to fetch the
   * output data produced by higher rank processes. When calling this
   * function, process 0 halts until a valid data chunk is transmitted
   * by the queried process.
   *
   * \param mpidata: `const mixMPI*` Pointer to a `mixMPI` instance.
   * \param pid: `int` Rank of the process that is transmitting data.
   * \return result: `int` An exit code (0 for success).
   */
  int mpireceive(const mixMPI* mpidata, int pid);

  /*! \brief Send output data to process 0.
   *
   * This function is invoked by non-zero ranked MPI processes when
   * they are ready to send back the output data. When a process meets
   * this function call, it halts until MPI process 0 asks for the
   * data transmission.
   *
   * \param mpidata: `const mixMPI*` Pointer to a `mixMPI` instance.
   * \param pid: `int` Rank of the process that is transmitting data.
   * \return result: `int` An exit code (0 for success).
   */
  int mpisend(const mixMPI *mpidata);
#endif // MPI_VERSION
};
// >>> END OF OUTPUT FOR INCLUSION <<<

// >>> OUTPUT FOR SPHERE <<<
/*! \brief Class to collect output information for scattering from a single sphere.
 *
 * The results of the calculation can be saved in different formats.
 * It is therefore convenient to have a proper memory structure that
 * allows for storing the results and flushing them in any of the
 * permitted formats with just one operation. The purpose of the
 * `SphereOutputInfo` class is to provide a wrapper for the output
 * of the spherical particle scattering solver.
 */
class SphereOutputInfo {
protected:
  //! \brief Flag for skipping mpisend() and mpireceive()
  int _skip_flag;
  //! \brief Number of incident azimuth calculations.
  int _num_theta;
  //! \brief Number of scattered azimuth calculations.
  int _num_thetas;
  //! \brief Number of incident elevation calculations.
  int _num_phi;
  //! \brief Number of scattered elevation calculations.
  int _num_phis;
  //! \brief ID of the first computed wavelength
  int _first_xi;
  
  /*! \brief Write the output to a HDF5 file.
   *
   * \param file_name: `const string &` Path to the output to be written.
   * \return result: `int` Exit code (0 if successful).
   */
  int write_hdf5(const std::string &file_name);
  
  /*! \brief Write the output to a legacy text file.
   *
   * This function takes care of writing the output using the legacy
   * formatted ASCII structure. If the output file does not exist, it
   * is created. If it exists, the new content is overwritten.
   *
   * \param output: `const string &` Path to the output to be written.
   * \return result: `int` Exit code (0 if successful).
   */
  int write_legacy(const std::string &output);
  
public:
  //! \brief Read-only view on skip_flag
  const int &skip_flag = _skip_flag;
  //! \brief Read-only view on the ID of the first scale
  const int &first_xi = _first_xi;
  //! \brief Number of spheres.
  int nsph;
  //! \brief Maximum field expansion order.
  int lm;
  //! \brief Incident polarization flag.
  int inpol;
  //! \brief Number of points for transition layer integration.
  int npnt;
  //! \brief Number of points for non-transition layer integration.
  int npntts;
  //! \brief Flag for reference to meridional plane.
  int isam;
  //! \brief Flag for dielectric function definition.
  int idfc;
  //! \brief First incident radiation azimuth angle.
  double th;
  //! \brief Incident radiation azimuth angle step.
  double thstp;
  //! \brief Last incident radiation azimuth angle.
  double thlst;
  //! \brief First scattered radiation azimuth angle.
  double ths;
  //! \brief Scattered radiation azimuth angle step.
  double thsstp;
  //! \brief Last scattered radiation azimuth angle.
  double thslst;
  //! \brief First incident radiation elevation angle.
  double ph;
  //! \brief Incident radiation elevation angle step.
  double phstp;
  //! \brief Last incident radiation elevation angle.
  double phlst;
  //! \brief First scattered radiation elevation angle.
  double phs;
  //! \brief Scattered radiation elevation angle step.
  double phsstp;
  //! \brief Last scattered radiation elevation angle.
  double phslst;
  //! \brief Number of directions to be explicitly solved.
  int ndirs;
  //! \brief Refractive index of external medium.
  double exri;
  //! \brief Number of scales (wavelengths)
  int nxi;
  //! \brief Number of scales handled by the current process.
  int xi_block_size;
  //! \brief Index of the wavelength for T-matrix output.
  int jwtm;
  //! \brief Number of sphere types.
  int configurations;
  //! \brief Highest expansion order achieved in calculations.
  int lcalc;
  //! \brief Harmonic functions argument.
  dcomplex arg;
  //! \brief Vector of scale (wavelength) indices.
  int *vec_jxi;
  //! \brief Vector of error severities (0 - success, 1 - DME).
  short *vec_ier;
  //! \brief Vector of vacuum wave numbers.
  double *vec_vk;
  //! \brief Vector of computed scales.
  double *vec_xi;
  //! \brief Vector of sphere sizes (one for every configuration and scale).
  double *vec_sphere_sizes;
  //! \brief Vector of sphere refractive indices (one for every configuration and scale).
  dcomplex *vec_sphere_ref_indices;
  //! \brief Vector of sphere scattering cross-sections.
  double *vec_scs;
  //! \brief Vector of sphere absorption cross-sections.
  double *vec_abs;
  //! \brief Vector of sphere extinction cross-sections.
  double *vec_exs;
  //! \brief Vector of sphere albedos.
  double *vec_albeds;
  //! \brief Vector of sphere scattering-to-geometric cross-sections.
  double *vec_scsrt;
  //! \brief Vector of sphere absorption-to-geometric cross-sections.
  double *vec_absrt;
  //! \brief Vector of sphere extinction-to-geometric cross-sections.
  double *vec_exsrt;
  //! \brief Vector of sphere forward scattering amplitudes.
  dcomplex *vec_fsas;
  //! \brief Vector of sphere QSCHU.
  double *vec_qschu;
  //! \brief Vector of sphere PSCHU.
  double *vec_pschu;
  //! \brief Vector of sphere S0MAG.
  double *vec_s0mag;
  //! \brief Vector of sphere average asymmetry parameter.
  double *vec_cosav;
  //! \brief Vector of sphere average radiation pressure force (N).
  double *vec_raprs;
  //! \brief Vector of sphere average extinction torque along incidence direction (parallel polarization).
  double *vec_tqek1;
  //! \brief Vector of sphere average extinction torque along incidence direction (perpendicular polarization).
  double *vec_tqek2;
  //! \brief Vector of sphere average scattering torque along incidence direction (parallel polarization).
  double *vec_tqsk1;
  //! \brief Vector of sphere average scattering torque along incidence direction (perpendicular polarization).
  double *vec_tqsk2;
  //! \brief Vector of total forward scattering amplitudes.
  dcomplex *vec_fsat;
  //! \brief Vector of total QSCHU.
  double *vec_qschut;
  //! \brief Vector of total PSCHU.
  double *vec_pschut;
  //! \brief Vector of total S0MAG.
  double *vec_s0magt;
  //! \brief Vector of incidence azimuth directions (one per incidence azimuth).
  double *vec_dir_tidg;
  //! \brief Vector of incidence elevation directions (one per incidence elevation).
  double *vec_dir_pidg;
  //! \brief Vector of scattering azimuth directions (one per scattering azimuth).
  double *vec_dir_tsdg;
  //! \brief Vector of scattering elevation directions (one per scattering elevation).
  double *vec_dir_psdg;
  //! \brief Vector of scattering angles (one per direction).
  double *vec_dir_scand;
  //! \brief Control parameter for incidence plane referred to meridional plane (one per direction).
  double *vec_dir_cfmp;
  //! \brief Control parameter for scattering plane referred to meridional plane (one per direction).
  double *vec_dir_sfmp;
  //! \brief Control parameter for incidence plane referred to scattering plane (one per direction).
  double *vec_dir_cfsp;
  //! \brief Control parameter for scattering plane referred to scattering plane (one per direction).
  double *vec_dir_sfsp;
  //! \brief Components of the unitary vector perpendicular to incidence plane (three per direction).
  double *vec_dir_un;
  //! \brief Components of the unitary vector perpendicular to scattering plane (three per direction).
  double *vec_dir_uns;
  //! \brief Vector of sphere differential scattering amplitude with polarization parallel to parallel incidence field.
  dcomplex *vec_dir_sas11;
  //! \brief Vector of sphere differential scattering amplitude with polarization perpendicular to the parallel incidence field.
  dcomplex *vec_dir_sas21;
  //! \brief Vector of sphere differential scattering amplitude with polarization perpendicular to perpendicular incidence field.
  dcomplex *vec_dir_sas12;
  //! \brief Vector of sphere differential scattering amplitude with polarization parallel the perpendicular incidence field.
  dcomplex *vec_dir_sas22;
  //! \brief Vector of differential radiation pressure force components along the X axis.
  double *vec_dir_fx;
  //! \brief Vector of differential radiation pressure force components along the Y axis.
  double *vec_dir_fy;
  //! \brief Vector of differential radiation pressure force components along the Z axis.
  double *vec_dir_fz;
  //! \brief Vector of sphere Mueller transormation matrices referred to meridional plane.
  double *vec_dir_muls;
  //! \brief Vector of sphere Mueller transormation matrices referred to scattering plane.
  double *vec_dir_mulslr;
  
  /*! \brief `SphereOutputInfo` default instance constructor.
   *
   * \param sc: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` instance.
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \param mpidata: `const mixMPI*` Pointer to a mixMPI instance.
   * \param first_xi: `int` Index of the first scale in output (optional, default is 1).
   * \param xi_length: `int` Number of scales tobe included in output (optional, default is 0, meaning all).
   */   
  SphereOutputInfo(
    ScattererConfiguration *sc, GeometryConfiguration *gc,
    const mixMPI *mpidata, int first_xi = 1, int xi_length = 0
  );

  /*! \brief `SphereOutputInfo` constructor from HDF5 input.
   *
   * \param hdf5_name: `const string &` Path to the HDF5 file to be read.
   */   
  SphereOutputInfo(const std::string &hdf5_name);

  /*! \brief `SphereOutputInfo` constructor for the dummy NULL case
   *
   * \param skip_flag: `const int` must be passed as the `1` constant.
   */   
  SphereOutputInfo(const int skip_flag);

  /*! \brief `InclusionOutputInfo` instance destroyer.
   */
  ~SphereOutputInfo();

  /*! \brief Estimate the size of the structure that would be built for given input.
   *
   * \param sc: `ScattererConfiguration *` Pointer to a `ScattererConfiguration` instance.
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \param first_xi: `int` Index of the first scale in output (optional, default is 1).
   * \param xi_length: `int` Number of scales tobe included in output (optional, default is all).
   * \return size: `long` Estimated instance size in bytes.
   */
  static long compute_size(
    ScattererConfiguration *sc, GeometryConfiguration *gc,
    int first_xi = 1, int xi_length = 0
  );
  
  /*! \brief Get the size of a `ClusterOutputInfo` instance in bytes.
   *
   * \return size: `long` Estimated instance size in bytes.
   */
  long compute_size();

  /*! \brief Insert in the current output data the data of another block.
   *
   * \param rhs: `const SphereOutputInfo &` Reference to the source data block.
   * \return result: `int` Exit code (0 if successful).
   */
  int insert(const SphereOutputInfo &rhs);

  /*! \brief Write the output to a file.
   *
   * \param output: `const string &` Path to the output to be written.
   * \param format: `const string &` Output format (one of LEGACY or HDF5).
   * \return result: `int` Exit code (0 if successful).
   */
  int write(const std::string &output, const std::string &format);

#ifdef MPI_VERSION
  /*! \brief Receive output data from worker processes.
   *
   * This function is invoked by the MPI rank-0 process to fetch the
   * output data produced by higher rank processes. When calling this
   * function, process 0 halts until a valid data chunk is transmitted
   * by the queried process.
   *
   * \param mpidata: `const mixMPI*` Pointer to a `mixMPI` instance.
   * \param pid: `int` Rank of the process that is transmitting data.
   * \return result: `int` An exit code (0 for success).
   */
  int mpireceive(const mixMPI *mpidata, int pid);

  /*! \brief Send output data to process 0.
   *
   * This function is invoked by non-zero ranked MPI processes when
   * they are ready to send back the output data. When a process meets
   * this function call, it halts until MPI process 0 asks for the
   * data transmission.
   *
   * \param mpidata: `const mixMPI*` Pointer to a `mixMPI` instance.
   * \param pid: `int` Rank of the process that is transmitting data.
   * \return result: `int` An exit code (0 for success).
   */
  int mpisend(const mixMPI *mpidata);
#endif // MPI_VERSION
};
// >>> END OF OUTPUT FOR SPHERE <<<

#endif // INCLUDE_OUTPUTS_H_
