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

/*! \file tfrfme.h
 *
 * \brief Representation of the trapping calculation objects.
 */

#ifndef INCLUDE_TFRFME_H_
#define INCLUDE_TFRFME_H_

/*! \brief Class to represent the first group of trapping swap data.
 */
class Swap1 {
protected:
  //! Index of the last element to be filled.
  int _last_index;
  //! Number of beam description wave-numbers.
  int _nkv;
  //! NLMMT = 2 * LM * (LM + 2)
  int _nlmmt;

  //! QUESTION: definition?
  dcomplex *_wk;

  /*! \brief Load a Swap1 instance from a HDF5 binary file.
   *
   * \param file_name: `string` Name of the file to be loaded.
   * \return instance: `Swap1 *` Pointer to a new Swap1 instance.
   */
  static Swap1 *from_hdf5(const std::string& file_name);

  /*! \brief Load a Swap1 instance from a legacy binary file.
   *
   * \param file_name: `string` Name of the file to be loaded.
   * \return instance: `Swap1 *` Pointer to a new Swap1 instance.
   */
  static Swap1 *from_legacy(const std::string& file_name);

  /*! \brief Save a Swap1 instance to a HDF5 binary file.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_hdf5(const std::string& file_name);

  /*! \brief Save a Swap1 instance to a legacy binary file.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_legacy(const std::string& file_name);

public:
  //! \brief Read only view on WK.
  const dcomplex *wk;
  
  /*! \brief Swap1 instance constructor.
   *
   * \param lm: `int` Maximum field expansion order.
   * \param nkv: `int` Number of beam description wave numbers.
   */
  Swap1(int lm, int nkv);

  /*! \brief Swap1 instance destroyer.
   */
  ~Swap1() { delete[] _wk; }

  /*! \brief Append an element at the end of the vector.
   *
   * \param value: `complex double` The value to be added to the vector.
   */
  void append(dcomplex value) { _wk[_last_index++] = value; }
  
  /*! \brief Load a Swap1 instance from binary file.
   *
   * \param file_name: `string` Name of the file.
   * \param mode: `string` Format of the file (can be either "HDF5"
   * or "LGEACY". Default is "LEGACY").
   * \return instance: `Swap1 *` Pointer to a newly created Swap1 instance.
   */
  static Swap1* from_binary(const std::string& file_name, const std::string& mode="LEGACY");

  /*! \brief Calculate the necessary amount of memory to create a new instance.
   *
   * \param lm: `int` Maximum field expansion order.
   * \param nkv: `int` Number of vector coordinates. QUESTION: correct?
   * \return size: `long` The necessary memory size in bytes.
   */
  static long get_size(int lm, int nkv);
  
  /*! \brief Bring the pointer to the next element at the start of vector.
   */
  void reset() { _last_index = 0; }

  /*! \brief Write a Swap1 instance to binary file.
   *
   * \param file_name: `string` Name of the file.
   * \param mode: `string` Format of the file (can be either "HDF5"
   * or "LGEACY". Default is "LEGACY").
   */
  void write_binary(const std::string& file_name, const std::string& mode="LEGACY");

  /*! \brief Test whether two instances of Swap1 are equal.
   *
   * \param other: `Swap1 &` Reference to the instance to be compared
   * with.
   * \return result: `bool` True, if the two instances are equal,
   * false otherwise.
   */
  bool operator ==(Swap1 &other);
};

/*! \brief Class to represent the second group of trapping swap data.
 */
class Swap2 {
protected:
  //! Index of the last vector element to be filled.
  int _last_vector;
  //! Index of the last matrix element to be filled.
  int _last_matrix;
  //! Number of beam description wave numbers.
  int _nkv;
  //! QUESTION: definition?
  double _apfafa;
  //! QUESTION: definition?
  double _pmf;
  //! QUESTION: definition?
  double _spd;
  //! QUESTION: definition?
  double _rir;
  //! QUESTION: definition?
  double _ftcn;
  //! QUESTION: definition?
  double _fshmx;
  //! QUESTION: definition?
  double _vxyzmx;
  //! Cartesian displacement. QUESTION: correct?
  double _delxyz;
  //! QUESTION: definition?
  double _vknmx;
  //! Wave number grid spacing.
  double _delk;
  //! Square of wave number grid spacing.
  double _delks;
  //! NLMMT = LM * (LM + 2) * 2
  int _nlmmt;
  //! Number of radial vector coordinates.
  int _nrvc;

  /*! \brief Load a Swap2 instance from a HDF5 binary file.
   *
   * \param file_name: `string` Name of the file to be loaded.
   * \return instance: `Swap2 *` Pointer to a new Swap2 instance.
   */
  static Swap2 *from_hdf5(const std::string& file_name);

  /*! \brief Load a Swap2 instance from a legacy binary file.
   *
   * \param file_name: `string` Name of the file to be loaded.
   * \return instance: `Swap2 *` Pointer to a new Swap2 instance.
   */
  static Swap2 *from_legacy(const std::string& file_name);

  /*! \brief Save a Swap2 instance to a HDF5 binary file.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_hdf5(const std::string& file_name);

  /*! \brief Save a Swap2 instance to a legacy binary file.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_legacy(const std::string& file_name);

public:
  //! Read-only view onthe index of the last vector element to be filled.
  const int &last_vector = _last_vector;
  //! Read-only view on the index of the last matrix element to be filled.
  const int &last_matrix = _last_matrix;
  //! Read-only view on the number of beam description wave numbers.
  const int &nkv = _nkv;
  //! QUESTION: definition?
  double *vkv;
  //! QUESTION: definition?
  double *vec_vkzm;
  //! QUESTION: definition?
  const double &apfafa = _apfafa;
  //! QUESTION: definition?
  const double &pmf = _pmf;
  //! QUESTION: definition?
  const double &spd = _spd;
  //! QUESTION: definition?
  const double &rir = _rir;
  //! QUESTION: definition?
  const double &ftcn = _ftcn;
  //! QUESTION: definition?
  const double &fshmx = _fshmx;
  //! QUESTION: definition?
  const double &vxyzmx = _vxyzmx;
  //! Cartesian displacement. QUESTION: correct?
  const double &delxyz = _delxyz;
  //! QUESTION: definition?
  const double &vknmx = _vknmx;
  //! QUESTION: definition?
  const double &delk = _delk;
  //! QUESTION: definition?
  const double &delks = _delks;
  //! NLMMT = LM * (LM + 2) * 2
  const int &nlmmt = _nlmmt;
  //! Read-only view on the number of radial vector coordinates.
  const int &nrvc = _nrvc;

  /*! \brief Swap2 instance constructor.
   *
   * \param nkv: `int` Number of beam description wave numbers.
   */
  Swap2(int nkv);

  /*! \brief Swap2 instance destroyer.
   */
  ~Swap2();
  
  /*! \brief Load a Swap2 instance from binary file.
   *
   * \param file_name: `string` Name of the file.
   * \param mode: `string` Format of the file (can be either "HDF5"
   * or "LGEACY". Default is "LEGACY").
   * \return instance: `Swap2 *` Pointer to a newly created Swap2 instance.
   */
  static Swap2* from_binary(const std::string& file_name, const std::string& mode="LEGACY");

  /*! \brief Calculate the necessary amount of memory to create a new instance.
   *
   * \param nkv: `int` Number of beam description wave numbers.
   * \return size: `long` The necessary memory size in bytes.
   */
  static long get_size(int nkv);

  /*! \brief Get the pointer to the VKV vector.
   *
   * \return value: `double *` Pointer to the VKV vector.
   */
  double *get_vector() { return vkv; }

  /*! \brief Append an element at the end of the matrix.
   *
   * \param value: `double` The value to be pushed in the matrix.
   */
  void push_matrix(double value);

  /*! \brief Append an element at the end of the vector.
   *
   * \param value: `double` The value to be pushed in the vector.
   */
  void push_vector(double value) { vkv[_last_vector++] = value; }

  /*! \brief Bring the matrix pointer to the start of the array.
   */
  void reset_matrix() { _last_matrix = 0; }

  /*! \brief Bring the vector pointer to the start of the array.
   */
  void reset_vector() { _last_vector = 0; }

  /*! \brief Set a parameter by its name and value.
   *
   * \param param_name: `string` Name of the parameter.
   * \param value: `double` The value of the parameter.
   */
  void set_param(const std::string& param_name, double value);

  /*! \brief Write a Swap2 instance to binary file.
   *
   * \param file_name: `string` Name of the file.
   * \param mode: `string` Format of the file (can be either "HDF5"
   * or "LGEACY". Default is "LEGACY").
   */
  void write_binary(const std::string& file_name, const std::string& mode="LEGACY");

  /*! \brief Test whether two instances of Swap2 are equal.
   *
   * \param other: `Swap1 &` Reference to the instance to be compared
   * with.
   * \return result: `bool` True, if the two instances are equal,
   * false otherwise.
   */
  bool operator ==(Swap2 &other);
};

/*! \brief Class to represent the trapping configuration.
 */
class TFRFME {
protected:
  //! NLMMT = 2 * LM * (LM + 2)
  int _nlmmt;
  //! NRVC = NXV * NYV * NZV
  int _nrvc;
  //! Beam description mode.
  int _lmode;
  //! Maximum field expansion order.
  int _lm;
  //! Number of beam description wave numbers.
  int _nkv;
  //! Number of computed X coordinates.
  int _nxv;
  //! Number of computed Y coordinates.
  int _nyv;
  //! Number of computed Z coordinates.
  int _nzv;
  //! Vacuum wave number.
  double _vk;
  //! External medium refractive index
  double _exri;
  //! Numerical aperture.
  double _an;
  //! Filling factor.
  double _ff;
  //! Lens transmission.
  double _tra;
  //! QUESTION: definition?
  double _spd;
  //! QUESTION: definition?
  double _frsh;
  //! QUESTION: definition?
  double _exril;
  //! Vector of computed x positions
  double *xv;
  //! Vector of computed y positions
  double *yv;
  //! Vector of computed z positions
  double *zv;

  /*! \brief Load a configuration instance from a HDF5 binary file.
   *
   * \param file_name: `string` Name of the file to be loaded.
   * \return instance: `TFRFME *` Pointer to a new trapping configuration
   * instance.
   */
  static TFRFME *from_hdf5(const std::string& file_name);

  /*! \brief Load a configuration instance from a legacy binary file.
   *
   * \param file_name: `string` Name of the file to be loaded.
   * \return instance: `TFRFME *` Pointer to a new trapping configuration
   * instance.
   */
  static TFRFME *from_legacy(const std::string& file_name);

  /*! \brief Save a configuration instance to a HDF5 binary file.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_hdf5(const std::string& file_name);

  /*! \brief Save a configuration instance to a legacy binary file.
   *
   * \param file_name: `string` Name of the file to be written.
   */
  void write_legacy(const std::string& file_name);

public:
  //! Read-only view on NLMMT.
  const int& nlmmt = _nlmmt;
  //! Read-only view on NRVC.
  const int& nrvc = _nrvc;
  //! Read-only view on field expansion mode identifier.
  const int& lmode = _lmode;
  //! Read-only view on maximum field expansion order.
  const int& lm = _lm;
  //! QUESTION: definition?
  const int& nkv = _nkv;
  //! Read-only view on number of computed X coordinates.
  const int& nxv = _nxv;
  //! Read-only view on number of computed Y coordinates.
  const int& nyv = _nyv;
  //! Read-only view on number of computed Z coordinates.
  const int& nzv = _nzv;
  //! Read-only view on vacuum wave number.
  const double& vk = _vk;
  //! Read-only view on external medium refractive index
  const double& exri = _exri;
  //! Read-only view on numeric aperture.
  const double& an = _an;
  //! Read-only view on filling factor.
  const double& ff = _ff;
  //! Read-only view on lens transmission.
  const double& tra = _tra;
  //! QUESTION: definition?
  const double& spd = _spd;
  //! QUESTION: definition?
  const double& frsh = _frsh;
  //! QUESTION: definition?
  const double& exril = _exril;
  //! QUESTION: definition?
  dcomplex *vec_wsum;
  
  /*! \brief Trapping configuration instance constructor.
   *
   * \param lmode: `int` Order expansion mode flag.
   * \param lm: `int` Maximum field expansion order.
   * \param nkv: `int` Number of wave vector coordinates. QUESTION: correct?
   * \param nxv: `int` Number of computed X coordinates.
   * \param nyv: `int` Number of computed Y coordinates.
   * \param nzv: `int` Number of computed Z coordinates.
   */
  TFRFME(int lmode, int lm, int nkv, int nxv, int nyv, int nzv);
  
  /*! \brief Trapping configuration instance destroyer.
   */
  ~TFRFME();

  /*! \brief Load a trapping configuration instance from binary file.
   *
   * \param file_name: `string` Name of the file.
   * \param mode: `string` Format of the file (can be either "HDF5"
   * or "LGEACY". Default is "LEGACY").
   * \return instance: `TFRFME *` Pointer to a newly created configuration
   * instance.
   */
  static TFRFME* from_binary(
    const std::string& file_name, const std::string& mode="LEGACY"
);

  /*! \brief Calculate the necessary amount of memory to create a new instance.
   *
   * \param lm: `int` Maximum field expansion order.
   * \param nkv: `int` Number of radial vector coordinates. QUESTION: correct?
   * \param nxv: `int` Number of computed X coordinates.
   * \param nyv: `int` Number of computed Y coordinates.
   * \param nzv: `int` Number of computed Z coordinates.
   * \return size: `long` The necessary memory size in bytes.
   */
  static long get_size(int lm, int nkv, int nxv, int nyv, int nzv);

  /*! \brief Get the pointer to the X coordinates vector.
   *
   * \return x: `double *` Pointer to X coordinates vector.
   */
  double *get_x() { return xv; }

  /*! \brief Get the pointer to the Y coordinates vector.
   *
   * \return y: `double *` Pointer to Y coordinates vector.
   */
  double *get_y() { return yv; }

  /*! \brief Get the pointer to the Z coordinates vector.
   *
   * \return z: `double *` Pointer to Z coordinates vector.
   */
  double *get_z() { return zv; }

  /*! \brief Set a configuration parameter.
   *
   * \param param_name: `string` Name of the parameter.
   * \param value: `double` Value to be stored as parameter.
   */
  void set_param(const std::string& param_name, double value);

  /*! \brief Write a trapping configuration instance to binary file.
   *
   * \param file_name: `string` Name of the file.
   * \param mode: `string` Format of the file (can be either "HDF5"
   * or "LGEACY". Default is "LEGACY").
   */
  void write_binary(const std::string& file_name, const std::string& mode="LEGACY");
  
  /*! \brief Test whether two instances of configuration are equal.
   *
   * \param other: `TFRFME &` Reference to the instance to be compared
   * with.
   * \return result: `bool` True, if the two instances are equal,
   * false otherwise.
   */
  bool operator ==(const TFRFME& other);
};
#endif
