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

/*! \file TransitionMatrix.h
 *
 * \brief Representation of the Transition Matrix.
 */

#ifndef INCLUDE_TRANSITIONMATRIX_H_
#define INCLUDE_TRANSITIONMATRIX_H_

/*! \brief Class to represent the Transition Matrix.
 */
class TransitionMatrix {
 protected:
  //! Matrix type identifier.
  int _is;
  //! Maximum field expansion order.
  int _l_max;
  //! Wave number in scale units.
  double _vk;
  //! External medium refractive index.
  double _exri;
  //! Sphere radius.
  double _sphere_radius;

  /*! \brief Build transition matrix from a HDF5 binary input file.
   *
   * This function takes care of the specific task of building a transition
   * matrix memory data structure from a HDF5 binary input file.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \return config: `TransitionMatrix *` Pointer to object containing the
   * transition matrix data.
   */
  static TransitionMatrix *from_hdf5(const std::string& file_name);
  
  /*! \brief Build transition matrix from a legacy binary input file.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \return config: `TransitionMatrix *` Pointer to object containing the
   * transition matrix data.
   */
  static TransitionMatrix *from_legacy(const std::string& file_name);

  /*! \brief Write the Transition Matrix to HDF5 binary output.
   *
   * This function takes care of the specific task of writing the transition
   * matrix memory data structure to a binary output file formatted according
   * to the HDF5 standard.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   */
  void write_hdf5(const std::string& file_name);
  
  /*! \brief Write transition matrix data to HDF5 binary output.
   *
   * This function takes care of the specific task of writing the transition
   * matrix memory data structure to a binary output file formatted according
   * to the HDF5 standard without a pre-existing instance. It is designed to
   * work for the case of a cluster of spheres.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \param nlemt: `np_int` Size of the matrix (2 * LE * (LE + 2)).
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param am0m: `complex double **`
   */
  static void write_hdf5(
    const std::string& file_name, np_int nlemt, int lm, double vk,
    double exri, dcomplex **am0m
  );
  
  /*! \brief Write transition matrix data to HDF5 binary output.
   *
   * This function takes care of the specific task of writing the transition
   * matrix memory data structure to a binary output file formatted according
   * to the HDF5 standard without a pre-existing instance. It is designed to
   * work for the case of a single sphere.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param rmi: `complex double **`
   * \param rei: `complex double **`
   * \param sphere_radius: `double` Radius of the sphere.
   */
  static void write_hdf5(
    const std::string& file_name, int lm, double vk, double exri,
    dcomplex **rmi, dcomplex **rei, double sphere_radius
  );

  /*! \brief Write the Transition Matrix to legacy binary output.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   */
  void write_legacy(const std::string& file_name);

  /*! \brief Write transition matrix data to binary output using legacy format.
   *
   * This function takes care of the specific task of writing the transition
   * matrix memory data structure to a binary output file formatted according
   * to the format used by the legacy FORTRAN code. It is designed to work for
   * the case of clusters of spheres.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \param nlemt: `np_int` Size of the matrix (2 * LE * (LE + 2)).
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param am0m: `complex double **`
   */
  static void write_legacy(
    const std::string& file_name, np_int nlemt, int lm, double vk,
    double exri, dcomplex **am0m
  );
  
  /*! \brief Write transition matrix data to binary output using legacy format.
   *
   * This function takes care of the specific task of writing the transition
   * matrix memory data structure to a binary output file formatted according
   * to the original code structure without a pre-existing instance. It is designed
   * to work for the case of a single sphere.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param rmi: `complex double **`
   * \param rei: `complex double **`
   * \param radius: `double` Radius of the sphere.
   */
  static void write_legacy(
    const std::string& file_name, int lm, double vk, double exri,
    dcomplex **rmi, dcomplex **rei, double radius
  );

 public:
  //! Read-only view on matrix type identifier.
  const int& is = _is;
  //! Read-only view on maximum field expansion order.
  const int& l_max = _l_max;
  //! Read-only view on wave number in scale units.
  const double& vk = _vk;
  //! Read-only view on external medium refractive index.
  const double& exri = _exri;
  //! Vectorized matrix elements.
  dcomplex *elements;
  //! Read-only view on sphere radius.
  const double& sphere_radius = _sphere_radius;
  //! Matrix shape
  int *shape;
  
  /*! \brief Default Transition Matrix instance constructor.
   *
   * \param is: `int` Matrix type identifier
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param _elems: `complex double *` Vectorized elements of the matrix.
   * \param radius: `double` Radius for the single sphere case (defaults to 0.0).
   */
  TransitionMatrix(
    int is, int lm, double vk, double exri, dcomplex *_elems, double radius=0.0
  );

  /*! \brief Transition Matrix instance constructor for single sphere.
   *
   * This constructor allocates the memory structure needed to represent the transition
   * matrix for the case of a single sphere.
   *
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param rmi: `complex double **`
   * \param rei: `complex double **`
   * \param radius: `double` Radius of the sphere.
   */
  TransitionMatrix(
     int lm, double vk, double exri, dcomplex **rmi, dcomplex **rei, double radius
  );

  /*! \brief Transition Matrix instance constructor for a cluster of spheres.
   *
   * This constructor allocates the memory structure needed to represent the transition
   * matrix for the case of a cluster of spheres.
   *
   * \param nlemt: `int` Size of the matrix (2 * LE * (LE + 2)).
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param am0m: `complex double **`
   */
  TransitionMatrix(int nlemt, int lm, double vk, double exri, dcomplex **am0m);

  /*! \brief Transition Matrix instance destroyer.
   */
  ~TransitionMatrix();
  
  /*! \brief Build transition matrix from binary input file.
   *
   * In some cases, it is necessary to perform calculations starting from a pre-computed
   * transition matrix. If this matrix is not available in memory (e.g. because it was
   * calculated by a different process), but it was saved in a binary file, this function
   * can be used to load it back in memory. The function operates in two modes: "LEGACY",
   * which reads the matrix data from a proprietary binary file, having the same structure
   * as the one used by the original FORTRAN code, and "HDF5", which, instead, reads the
   * data from an input file complying with the HDF5 format.
   *
   * \param file_name: `string` Name of the binary configuration data file.
   * \param mode: `string` Binary encoding. Can be one of ["LEGACY", "HDF5"]. Optional
   * (default is "LEGACY").
   * \return config: `TransitionMatrix *` Pointer to object containing the transition
   * matrix data.
   */
  static TransitionMatrix* from_binary(const std::string& file_name, const std::string& mode="LEGACY");

  /*! \brief Write the Transition Matrix to a binary file.
   *
   * This function writes a hard-copy of the transition matrix to an output file, making
   * it available for subsequent processes to reload. The function operates in two modes:
   * "LEGACY", which writes a proprietary binary file, using the same structure of the
   * original FORTRAN code, and "HDF5", which, instead, writes the output to a file using
   * the HDF5 format, thus leaving it available for inspection with external tools.
   *
   * \param file_name: `string` Name of the file to be written.
   * \param mode: `string` Binary encoding. Can be one of ["LEGACY", "HDF5"] . Optional
   * (default is "LEGACY").
   */
  void write_binary(const std::string& file_name, const std::string& mode="LEGACY");
  
  /*! \brief Write a cluster Transition Matrix to a binary file without instanciating it.
   *
   * Transition Matrix data can take a large amount of memory. For such reason, attempts
   * to create TransitionMatrix instances only for writing purposes can create
   * unnecessary resource consumption and computing time to duplicate the data into
   * the output buffer. This function offers output to file as a static method. It
   * takes the arguments of a constructor together with the usual arguments to specify
   * the output file name and format, to write the required data directly to a file,
   * without creating a new TransitionMatrix instance. The implementation works for
   * TransitionMatrix objects built for the CLUSTER case. It belongs to the public class
   * interface and it calls the proper versions of `write_legacy()` and `write_hdf5()`,
   * depending on the requested output format.
   * 
   * \param file_name: `string` Name of the file to be written.
   * \param nlemt: `np_int` Size of the matrix (2 * LE * (LE + 2)).
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param am0m: `complex double **`
   * \param mode: `string` Binary encoding. Can be one of ["LEGACY", "HDF5"] . Optional
   * (default is "LEGACY").
   */
  static void write_binary(
     const std::string& file_name, np_int nlemt, int lm, double vk,
     double exri, dcomplex **am0m, const std::string& mode="LEGACY"
  );
  
  /*! \brief Write a single sphere Transition Matrix to a binary file without instanciating it.
   *
   * Transition Matrix data can take a large amount of memory. For such reason, attempts
   * to create TransitionMatrix instances only for writing purposes can create
   * unnecessary resource consumption and computing time to duplicate the data into
   * the output buffer. This function offers output to file as a static method. It
   * takes the arguments of a constructor together with the usual arguments to specify
   * the output file name and format, to write the required data directly to a file,
   * without creating a new TransitionMatrix instance. The implementation works for
   * TransitionMatrix objects built for the single sphere case. It belongs to the public
   * class interface and it calls the proper versions of `write_legacy()` and `write_hdf5()`,
   * depending on the requested output format.
   *
   * \param file_name: `string` Name of the file to be written.
   * \param lm: `int` Maximum field expansion order.
   * \param vk: `double` Wave number in scale units.
   * \param exri: `double` External medium refractive index.
   * \param rmi: `complex double **`
   * \param rei: `complex double **`
   * \param sphere_radius: `double` Radius of the sphere.
   * \param mode: `string` Binary encoding. Can be one of ["LEGACY", "HDF5"] . Optional
   * (default is "LEGACY").
   */
  static void write_binary(
    const std::string& file_name, int lm, double vk, double exri,
    dcomplex **rmi, dcomplex **rei, double sphere_radius,
    const std::string& mode="LEGACY"
  );

  /*! \brief Test whether two instances of TransitionMatrix are equal.
   *
   * Transition matrices can be the result of a calculation or of a file input operation,
   * reading from a previously computed object. The `==` operator tests for the equivalence
   * of two transition matrices, returning `true` if they are equivalnet, `false` otherwise.
   *
   * \param other: `TransitionMatrix &` Reference to the instance to be compared with.
   * \return result: `bool` True, if the two instances are equal, false otherwise.
   */
  bool operator ==(TransitionMatrix &other);
};

#endif
