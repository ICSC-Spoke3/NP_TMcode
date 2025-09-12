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

/*! \file Configuration.h
 *
 * \brief Configuration data structures.
 *
 * To handle the calculation of a scattering problem, the code needs a set
 * of configuration parameters that define the properties of the scattering
 * particle, of the incident radiation field and of the geometry of the
 * problem. The necessary information is managed through the use of two
 * classes: `ScattererConfiguration` and `GeometryConfiguration`. The first
 * class, `ScattererConfiguration`, undertakes the role of describing the
 * scattering particle properties, by defining its structure in terms of a
 * distribution of spherical sub-particles, and the optical properties of
 * the constituting materials. It also defines the scaling vector, which
 * tracks the ratio of particle size and radiation wavelength, eventually
 * allowing for the iteration of the calculation on different wavelengths.
 * Multiple materials and layered structures are allowed, while sphere
 * compenetration is not handled, due to the implied discontinuity on the
 * spherical symmetry of the elementary sub-particles. The second class,
 * `GeometryConfiguration`, on the contrary, describes the properties of
 * the radiation field impinging on the particle. It defines the incoming
 * radiation polarization state, the incident direction and the scattering
 * directions that need to be accounted for. Both classes expose methods to
 * read the required configuration data from input files formatted according
 * to the same rules expected by the original code. In addition, they offer
 * perform I/O operations towards proprietary and standard binary formats.
 */

#ifndef INCLUDE_CONFIGURATION_H_
#define INCLUDE_CONFIGURATION_H_

class mixMPI;

/**
 * \brief A class to represent the configuration of the scattering geometry.
 *
 * GeometryConfiguration is a class designed to store the necessary configuration
 * data to describe the scattering geometry, including the distribution of the
 * particle components, the orientation of the incident and scattered radiation
 * fields and their polarization properties.
 */
class GeometryConfiguration {

protected:
  //! \brief Number of spherical components.
  int _number_of_spheres;
  //! \brief Maximum field expansion order.
  int _l_max;
  //! \brief Maximum internal field expansion order.
  int _li;
  //! \brief Maximum external field expansion order.
  int _le;
  //! \brief Maximum dimension of allocated matrix allowance (deprecated).
  np_int _mxndm;
  //! \brief Flag for intensity.
  int _iavm;
  //! \brief Incident field polarization status (0 - linear, 1 - circular).
  int _in_pol;
  //! \brief Number of points for transition layer integration.
  int _npnt;
  //! \brief Number of points for non-transition layer integration.
  int _npntts;
  //! \brief Type of meridional plane definition.
  int _isam;
  //! \brief Scale index of the T-matrix output.
  int _jwtm;
  //! \brief Incident field initial azimuth.
  double _in_theta_start;
  //! \brief Incident field azimuth step.
  double _in_theta_step;
  //! \brief Incident field final azimuth.
  double _in_theta_end;
  //! \brief Scattered field initial azimuth.
  double _sc_theta_start;
  //! \brief Scattered field azimuth step.
  double _sc_theta_step;
  //! \brief Scattered field final azimuth.
  double _sc_theta_end;
  //! \brief Incident field initial elevation.
  double _in_phi_start;
  //! \brief Incident field elevation step.
  double _in_phi_step;
  //! \brief Incident field final elevation.
  double _in_phi_end;
  //! \brief Scattered field initial elevation.
  double _sc_phi_start;
  //! \brief Scattered field elevation step.
  double _sc_phi_step;
  //! \brief Scattered field final elevation.
  double _sc_phi_end;
  //! \brief Vector of spherical components X coordinates.
  double *_sph_x;
  //! \brief Vector of spherical components Y coordinates.
  double *_sph_y;
  //! \brief Vector of spherical components Z coordinates.
  double *_sph_z;
  //! \brief Flag for matrix inversion refinement.
  short _refine_flag;
  //! \brief Flag for dynamic order management.
  short _dyn_order_flag;

public:
  //! \brief Read-only view on number of spherical components.
  const int& number_of_spheres = _number_of_spheres;
  //! \brief Read-only view on maximum field expansion order.
  const int& l_max = _l_max;
  //! \brief Read-only view on maximum internal field expansion order.
  const int& li = _li;
  //! \brief Read-only view on maximum external field expansion order.
  const int& le = _le;
  //! \brief Read-only view on maximum dimension of allocated matrix allowance (deprecated).
  const np_int& mxndm = _mxndm;
  //! \brief Read-only view on the intensity mode flag.
  const int& iavm = _iavm;
  //! \brief Read-only view on incident field polarization status (0 - linear, 1 - circular).
  const int& in_pol = _in_pol;
  //! \brief Read-only view on number of points for transition layer integration.
  const int& npnt = _npnt;
  //! \brief Read-only view on number of points for non-transition layer integration.
  const int& npntts = _npntts;
  //! \brief Read-only view on type of meridional plane definition.
  const int& isam = _isam;
  //! \brief Read-only view on scale index for T-matrix output.
  const int& jwtm = _jwtm;
  //! \brief Read-only view on incident field initial azimuth.
  const double& in_theta_start = _in_theta_start;
  //! \brief Read-only view on incident field azimuth step.
  const double& in_theta_step = _in_theta_step;
  //! \brief Read-only view on incident field final azimuth.
  const double& in_theta_end = _in_theta_end;
  //! \brief Read-only view on scattered field initial azimuth.
  const double& sc_theta_start = _sc_theta_start;
  //! \brief Read-only view on scattered field azimuth step.
  const double& sc_theta_step = _sc_theta_step;
  //! \brief Read-only view on scattered field final azimuth.
  const double& sc_theta_end = _sc_theta_end;
  //! \brief Read-only view on incident field initial elevation.
  const double& in_phi_start = _in_phi_start;
  //! \brief Read-only view on incident field elevation step.
  const double& in_phi_step = _in_phi_step;
  //! \brief Read-only view on incident field final elevation.
  const double& in_phi_end = _in_phi_end;
  //! \brief Read-only view on scattered field initial elevation.
  const double& sc_phi_start = _sc_phi_start;
  //! \brief Read-only view on scattered field elevation step.
  const double& sc_phi_step = _sc_phi_step;
  //! \brief Read-only view on scattered field final elevation.
  const double& sc_phi_end = _sc_phi_end;
  //! \brief Read-only view on flag for matrix inversion refinement.
  const short& refine_flag = _refine_flag;
  //! \brief Read-only view on flag for dynamic order management.
  const short& dyn_order_flag = _dyn_order_flag;
  
  /*! \brief Build a scattering geometry configuration structure.
   *
   * \param nsph: `int` Number of spheres to be used in calculation.
   * \param lm: `int` Maximum field angular momentum expansion order.
   * \param in_pol: `int` Incident field polarization status
   * \param npnt: `int` Number of transition points. QUESTION: correct?
   * \param npntts: `int` Transition smoothness. QUESTION: correct?
   * \param meridional_type: `int` Type of meridional plane definition (<0
   * for incident angles, 0 if determined by incidence and observation, =1
   * accross z-axis for incidence and observation, >1 across z-axis as a
   * function of incidence angles for fixed scattering).
   * \param li: `int`
   * \param le: `int`
   * \param mxndm: `int`
   * \param iavm: `int`
   * \param x: `double*` Vector of spherical components X coordinates.
   * \param y: `double*` Vector of spherical components Y coordinates.
   * \param z: `double*` Vector of spherical components Z coordinates.
   * \param in_th_start: `double` Incident field starting azimuth angle.
   * \param in_th_step: `double` Incident field azimuth angle step.
   * \param in_th_end: `double` Incident field final azimuth angle.
   * \param sc_th_start: `double` Scattered field starting azimuth angle.
   * \param sc_th_step: `double` Scattered field azimuth angle step.
   * \param sc_th_end: `double` Scattered field final azimuth angle.
   * \param in_ph_start: `double` Incident field starting elevation angle.
   * \param in_ph_step: `double` Incident field elevation angle step.
   * \param in_ph_end: `double` Incident field final elevation angle.
   * \param sc_ph_start: `double` Scattered field starting elevation angle.
   * \param sc_ph_step: `double` Scattered field elevation angle step.
   * \param sc_ph_end: `double` Scattered field final elevation angle.
   * \param jwtm: `int` Transition Matrix layer ID. QUESTION: correct?
   */
  GeometryConfiguration(
			int nsph, int lm, int in_pol, int npnt, int npntts, int meridional_type,
			int li, int le, np_int mxndm, int iavm,
			double *x, double *y, double *z,
			double in_th_start, double in_th_step, double in_th_end,
			double sc_th_start, double sc_th_step, double sc_th_end,
			double in_ph_start, double in_ph_step, double in_ph_end,
			double sc_ph_start, double sc_ph_step, double sc_ph_end,
			int jwtm
			);
  
  /*! \brief Build a scattering geometry configuration structure copying it from an existing one.
   *
   * \param rhs: `GeometryConfiguration` preexisting object to copy from.
   */
  GeometryConfiguration(const GeometryConfiguration& rhs);

#ifdef MPI_VERSION
  /*! \brief Build a scattering geometry configuration structure copying it via MPI from MPI process 0.
   *
   * \param rhs: `mixMPI *` pointer to the mpidata instance to use for the MPI communications.
   */
  GeometryConfiguration(const mixMPI *mpidata);

  /*! \brief Broadcast over MPI the GeometryConfiguration instance from MPI process 0 to all others.
   *
   * When using MPI, the initial GeometryConfiguration instance created by MPI process 0
   * needs to be replicated on all other processes. This function sends it using
   * MPI broadcast calls. The MPI broadcast calls in this function must match those
   * in the constructor using the mixMPI pointer.
   *
   * \param mpidata: `mixMPI *` Pointer to the mpi structure used to do the MPI broadcast.
   */
  void mpibcast(const mixMPI *mpidata);
#endif

  /*! \brief Destroy a GeometryConfiguration instance.
   */
  ~GeometryConfiguration();

  /*! \brief Build geometry configuration from legacy configuration input file.
   *
   * To allow for consistency tests and backward compatibility, geometry
   * configurations can be built from legacy configuration files. This function
   * replicates the approach implemented by the FORTRAN SPH and CLU codes, but
   * using a C++ oriented work-flow.
   *
   * \param file_name: `string` Name of the legacy configuration data file.
   * \return config: `GeometryConfiguration*` Pointer to object containing the
   * configuration data.
   */
  static GeometryConfiguration *from_legacy(const std::string& file_name);
  
  /*! \brief Get the X coordinate of a sphere by its index.
   *
   * This is a specialized function to access the X coordinate of a sphere through
   * its index.
   *
   * \param index: `int` Index of the scale to be retrieved.
   * \return scale: `double` The X coordinate of the requested sphere.
   */
  double get_sph_x(int index) { return _sph_x[index]; }
  
  /*! \brief Get the Y coordinate of a sphere by its index.
   *
   * This is a specialized function to access the Y coordinate of a sphere through
   * its index.
   *
   * \param index: `int` Index of the scale to be retrieved.
   * \return scale: `double` The Y coordinate of the requested sphere.
   */
  double get_sph_y(int index) { return _sph_y[index]; }
  
  /*! \brief Get the Z coordinate of a sphere by its index.
   *
   * This is a specialized function to access the Z coordinate of a sphere through
   * its index.
   *
   * \param index: `int` Index of the scale to be retrieved.
   * \return scale: `double` The Z coordinate of the requested sphere.
   */
  double get_sph_z(int index) { return _sph_z[index]; }

};

/**
 * \brief A class to represent scatterer configuration objects.
 *
 * ScattererConfiguration is a class designed to store the necessary configuration
 * data to describe the scatterer properties.
 */
class ScattererConfiguration {

protected:
  //! \brief Matrix of dielectric parameters with size [NON_TRANS_LAYERS x N_SPHERES x N_SCALES].
  dcomplex ***_dc0_matrix;
  //! \brief Vector of sphere radii expressed in m, with size [N_SPHERES].
  double *_radii_of_spheres;
  //! \brief Vector of sphere ID numbers, with size [N_SPHERES].
  int *_iog_vec;
  //! \brief Vector of layer numbers for every sphere, with size [CONFIGURATIONS].
  int *_nshl_vec;
  //! \brief Vector of scale parameters, with size [N_SCALES].
  double *_scale_vec;
  //! \brief Name of the reference variable type (one of XIV, WNS, WLS, PUS, EVS).
  std::string _reference_variable_name;
  //! \brief Number of spherical components.
  int _number_of_spheres;
  //! \brief Number of configurations.
  int _configurations;
  //! \brief Number of scales to use in calculation.
  int _number_of_scales;
  //! \brief Type of dielectric functions (<0 at XIP, =0 as function of XI, >0 constants).
  int _idfc;
  //! \brief External medium dielectric constant. QUESTION: correct?
  double _exdc;
  //! \brief WP. QUESTION: better definition?
  double _wp;
  //! \brief Peak XI. QUESTION: correct?
  double _xip;
  //! \brief Maximum number of layers for the particle components.
  int _max_layers;
  //! \brief Flag to control whether to add an external layer.
  bool _use_external_sphere;

  /*! \brief Build configuration from a HDF5 binary input file.
   *
   * This is the function called by the public method `from_binary()` in case of
   * HDF5 mode selection. This function creates a configuration structure from
   * a binary file written according to the HDF5 format standard.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \return config: `ScattererConfiguration*` Pointer to object containing the
   * scatterer configuration data.
   */
  static ScattererConfiguration *from_hdf5(const std::string& file_name);
  
  /*! \brief Build configuration from legacy binary input file.
   *
   * This is the function called by the public method `from_binary()` in case of
   * legacy mode selection. This function creates a configuration structure from
   * a binary file written according to the proprietary mode used by the original
   * FORTRAN code.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \return config: `ScattererConfiguration*` Pointer to object containing the
   * scatterer configuration data.
   */
  static ScattererConfiguration *from_legacy(const std::string& file_name);

  /*! \brief Write the scatterer configuration data to HDF5 binary output.
   *
   * This function is invoked by the public method `write_binary()` with the
   * "HDF5" format mode. It undertakes the task of writing the configuration
   * information to a binary file using the standard HDF5 format.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   */
  void write_hdf5(const std::string& file_name);
  
  /*! \brief Write the scatterer configuration data to legacy binary output.
   *
   * This function is invoked by the public method `write_binary()` with the
   * "LEGACY" format mode. It undertakes the task of writing the configuration
   * information to a binary file using a proprietary format, as it was done
   * originally in the FORTRAN code.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   */
  void write_legacy(const std::string& file_name);
public:
  //! \brief Read-only view on name of the reference variable type.
  const std::string& reference_variable_name = _reference_variable_name;
  //! \brief Read-only view on number of spherical components.
  const int& number_of_spheres = _number_of_spheres;
  //! \brief Read-only view on number of configurations.
  const int& configurations = _configurations;
  //! \brief Read-only view on number of scales to use in calculation.
  const int& number_of_scales = _number_of_scales;
  //! \brief Read-only view on type of dielectric functions.
  const int& idfc = _idfc;
  //! \brief Read-only view on external medium dielectric constant.
  const double& exdc = _exdc;
  //! \brief Read-only view on WP.
  const double& wp = _wp;
  //! \brief Read-only view on peak XI.
  const double& xip = _xip;
  //! \brief Read-only view on the maximum number of layers for the particle components.
  const int& max_layers = _max_layers;
  //! \brief Read-only view on flag to control whether to add an external layer.
  const bool& use_external_sphere = _use_external_sphere;
  //! \brief Matrix of fractional transition radii with size [CONFIGURATIONS x LAYERS].
  double **_rcf;
  
  /*! \brief Build a scatterer configuration structure.
   *
   * Prepare a default configuration structure by allocating the necessary
   * memory structures.
   *
   * \param nsph: `int` The number of spheres in the simulation.
   * \param configs: `int` Number of spherical monometer configuration types.
   * \param scale_vector: `double*` The radiation-particle scale vector.
   * \param nxi: `int` The number of radiation-particle scalings.
   * \param variable_name: `string` The name of the radiation-particle scaling type.
   * \param iog_vector: `int*` Array of sphere identification numbers. QUESTION: correct?
   * \param ros_vector: `double*` Sphere radius array.
   * \param nshl_vector: `int*` Array of layer numbers.
   * \param rcf_vector: `double**` Array of fractional break radii. QUESTION: correct?
   * \param dielectric_func_type: `int` Type of dielectric function definition (=0 for constant,
   * \>0 as function of scale parameter, <0 for functions at XIP value and XI is scale factor
   * for dimensions).
   * \param dc_matrix: `complex double ***` Matrix of reference dielectric constants.
   * \param has_external: `bool` Flag to set whether to add an external spherical layer.
   * \param exdc: `double` External medium dielectric constant.
   * \param wp: `double` wp
   * \param xip: `double` xip
   */
  ScattererConfiguration(
			 int nsph,
			 int configs,
			 double *scale_vector,
			 int nxi,
			 const std::string& variable_name,
			 int *iog_vector,
			 double *ros_vector,
			 int *nshl_vector,
			 double **rcf_vector,
			 int dielectric_func_type,
			 dcomplex ***dc_matrix,
			 bool has_external,
			 double exdc,
			 double wp,
			 double xip
  );
  
  /*! \brief Build a scatterer configuration structure copying its contents from a preexisting one.
   *
   * Prepare a default configuration structure by allocating the necessary
   * memory structures.
   *
   * \param rhs: `ScattererConfiguration&` Reference to the ScattereConfiguration
   * object to be copied.
   */
  ScattererConfiguration(const ScattererConfiguration& rhs);

#ifdef MPI_VERSION
  /*! \brief Build a scatterer configuration structure copying it via MPI from MPI process 0.
   *
   * \param rhs: `mixMPI *` pointer to the mpidata instance to use for the MPI communications.
   */
  ScattererConfiguration(const mixMPI *mpidata);

  /*! \brief Broadcast over MPI the ScattererConfiguration instance from MPI process 0 to all others.
   *
   * When using MPI, the initial ScattererConfiguration instance created by MPI process 0
   * needs to be replicated on all other processes. This function sends it using
   * MPI broadcast calls. The MPI broadcast calls in this function must match those
   * in the constructor using the mixMPI pointer.
   *
   * \param mpidata: `mixMPI *` Pointer to the mpi structure used to do the MPI broadcast.
   */
  void mpibcast(const mixMPI *mpidata);
#endif

  /*! \brief Destroy a scatterer configuration instance.
   */
  ~ScattererConfiguration();

  /*! \brief Build configuration from binary configuration input file.
   *
   * The configuration step can save configuration data as a binary file. The original
   * FORTRAN code used this possibility to manage communication between the configuring
   * code and the calculation program. This possibility is maintained, in case the
   * configuration step needs to be separated from the calculation execution. In this
   * case, `from_binary()` is the class method that restores a ScattererConfiguration
   * object from a previously saved binary file.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \param mode: `string` Binary encoding. Can be one of "LEGACY", ... . Optional
   * (default is "LEGACY").
   * \return config: `ScattererConfiguration*` Pointer to object containing the
   * scatterer configuration data.
   */
  static ScattererConfiguration* from_binary(const std::string& file_name, const std::string& mode = "LEGACY");

  /*! \brief Build scatterer configuration from legacy configuration input file.
   *
   * To allow for consistency tests and backward compatibility, ScattererConfiguration
   * objects can be built from legacy configuration files. This function replicates
   * the approach implemented by the FORTRAN EDFB code, but using a C++ oriented
   * work-flow.
   *
   * \param file_name: `string` Name of the legacy configuration data file.
   * \return config: `ScattererConfiguration*` Pointer to object containing the
   * scatterer configuration data.
   */
  static ScattererConfiguration* from_dedfb(const std::string& file_name);

  /*! \brief Get the dielectric constant of a material for a specific wavelength.
   *
   * Dielectric constants are stored in a 3D complex matrix, whose dimensions map
   * to [NUMBER_OF_CONFIGURATIONS x NUMBER_OF_SPHERES x NUMBER_OF_SCALES]. This
   * function extracts such values from the matrix through their indices.
   *
   * \param i: `int` Index of the configuration.
   * \param j: `int` Index of the sphere.
   * \param k: `int` Index of the current scale.
   * \return radius: `dcomplex` The requested dielectric constant.
   */
  dcomplex get_dielectric_constant(int i, int j, int k) { return _dc0_matrix[i][j][k]; }
  
  /*! \brief Get the ID of a configuration from the index of the sphere.
   *
   * This is a specialized function to get a configuration ID through the index of
   * the sphere it applies to.
   *
   * \param index: `int` Index of the sphere.
   * \return ID: `int` ID of the configuration to be applied.
   */
  int get_iog(int index) { return _iog_vec[index]; }
  
  /*! \brief Get the maximum radius of the sphere components.
   *
   * \return radius: `double` The radius of the largest sphere.
   */
  double get_max_radius();
  
  /*! \brief Get the number of layers for a given configuration.
   *
   * This is a specialized function to get the number of layers in a specific
   * configuration.
   *
   * \param index: `int` Index of the configuration.
   * \return nl: `int` The number of layers for the given configuration.
   */
  int get_nshl(int index) { return _nshl_vec[index]; }

  /*! \brief Get the radius of the smallest sphere containing the particle.
   *
   * \param gc: `GeometryConfiguration *` Pointer to a `GeometryConfiguration` instance.
   * \return radius: `double` The radius of the sphere containing the particle.
   */
  double get_particle_radius(GeometryConfiguration *gc);
  
  /*! \brief Get the radius of a sphere by its index.
   *
   * This is a specialized function to get the radius of a sphere through its
   * index.
   *
   * \param index: `int` Index of the ID to be retrieved.
   * \return radius: `double` The requested sphere radius.
   */
  double get_radius(int index) { return _radii_of_spheres[index]; }
  
  /*! \brief Get the value of a scale by its index.
   *
   * This is a specialized function to access a scale (generally a wavelength),
   * through its index.
   *
   * \param row: `int` Row index of the element to be retrieved.
   * \param column: `int` Column index of the element to be retrieved.
   * \return scale: `double` The desired scale.
   */
  double get_rcf(int row, int column) { return _rcf[row][column]; }

  /*! \brief Get the value of a scale by its index.
   *
   * This is a specialized function to access a scale (generally a wavelength),
   * through its index.
   *
   * \param index: `int` Index of the scale to be retrieved.
   * \return scale: `double` The desired scale.
   */
  double get_scale(int index) { return _scale_vec[index]; }
  
  /*! \brief Print the contents of the configuration object to terminal.
   *
   * In case of quick debug testing, `ScattererConfiguration.print()` allows printing
   * a formatted summary of the configuration data to terminal.
   */
  void print();

  /*! \brief Write the scatterer configuration data to binary output.
   *
   * The execution work-flow may be split in a configuration step and one or more
   * calculation steps. In case the calculation is not being run all-in-one, it can
   * be useful to save the configuration data. `ScattererConfiguration.write_binary()`
   * performs the operation of saving the configuration in binary format. This function
   * can work in legacy mode, to write backward compatible configuration files, as well
   * as by wrapping the data into common scientific formats.
   *
   * \param file_name: `string` Name of the file to be written.
   * \param mode: `string` Binary encoding. Can be one of ["LEGACY", "HDF5"] . Optional
   * (default is "LEGACY").
   */
  void write_binary(const std::string& file_name, const std::string& mode="LEGACY");

  /*! \brief Write the scatterer configuration data to formatted text output.
   *
   * Writing configuration to formatted text is an optional operation, which may turn
   * out to be useful for consistency checks. As a matter of fact, formatted configuration
   * output is not read back by the FORTRAN code work-flow and it can be safely omitted,
   * unless there is a specific interest in assessing that the legacy code and the
   * updated one are doing the same thing.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_formatted(const std::string& file_name);

  /*! \brief Test whether two instances of ScattererConfiguration are equal.
   *
   * ScattererConfiguration objects can be obtained in a variety of manner. They can
   * be constructed manually, through the class constructor, they can be read from
   * formatted configuration files, or they can be loaded from binary files. The `==`
   * operator tests for the equivalence of two configuration instances, returning `true`
   * if they are equivalnet, `false` otherwise.
   *
   * \param other: `ScattererConfiguration &` Reference to the instance to be compared with.
   * \return result: `bool` True, if the two instances are equal, false otherwise.
   */
  bool operator ==(const ScattererConfiguration &other);

};

#endif
