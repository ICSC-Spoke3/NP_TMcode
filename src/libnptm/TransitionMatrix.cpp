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

/*! \file TransitionMatrix.cpp
 *
 * \brief Implementation of the Transition Matrix structure.
 * \ingroup model
 */
#include <exception>
#include <fstream>
#include <hdf5.h>
#include <string>

#ifdef USE_MPI
#ifndef MPI_VERSION
#include <mpi.h>
#endif
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

#ifndef INCLUDE_TRANSITIONMATRIX_H_
#include "../include/TransitionMatrix.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_UTILS_H_
#include "../include/utils.h"
#endif

using namespace std;


TransitionMatrix::~TransitionMatrix() {
  if (elements != NULL) delete[] elements;
  if (shape != NULL) delete[] shape;
}

TransitionMatrix::TransitionMatrix(
  int is, int lm, double vk, double exri, dcomplex *elems, double radius
) {
  _is = is;
  _l_max = lm;
  _vk = vk;
  _exri = exri;
  elements = elems;
  _sphere_radius = radius;
  shape = new int[2]();
  if (is == 1111) {
    shape[0] = l_max;
    shape[1] = 2;
  } else if (is == 1) {
    const int nlemt = 2 * l_max * (l_max + 2);
    shape[0] = nlemt;
    shape[1] = nlemt;
  }
}

TransitionMatrix::TransitionMatrix(
  int lm, double vk, double exri, dcomplex **rmi,
  dcomplex **rei, double radius
) {
  _is = 1111;
  shape = new int[2];
  shape[0] = lm;
  shape[1] = 2;
  _l_max = lm;
  _vk = vk;
  _exri = exri;
  _sphere_radius = radius;
  elements = new dcomplex[2 * _l_max]();
  for (int ei = 0; ei < _l_max; ei++) {
    elements[2 * ei] = -1.0 / rmi[ei][0];
    elements[2 * ei + 1] = -1.0 / rei[ei][0];
  }
}

TransitionMatrix::TransitionMatrix(
  int nlemt, int lm, double vk, double exri, dcomplex **am0m
) {
  _is = 1;
  shape = new int[2];
  shape[0] = nlemt;
  shape[1] = nlemt;
  _l_max = lm;
  _vk = vk;
  _exri = exri;
  _sphere_radius = 0.0;
  elements = new dcomplex[nlemt * nlemt]();
  for (int ei = 0; ei < nlemt; ei++) {
    for (int ej = 0; ej < nlemt; ej++) elements[nlemt * ei + ej] = am0m[ei][ej];
  }
}

TransitionMatrix* TransitionMatrix::from_binary(const std::string& file_name, const std::string& mode) {
  TransitionMatrix *tm = NULL;
  if (mode.compare("LEGACY") == 0) {
    tm = TransitionMatrix::from_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    tm = TransitionMatrix::from_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"";
    throw UnrecognizedFormatException(message);
  }
  return tm;
}

TransitionMatrix* TransitionMatrix::from_hdf5(const std::string& file_name) {
  TransitionMatrix *tm = NULL;
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(file_name, flags);
  herr_t status = hdf_file->get_status();
  if (status == 0) {
    int is;
    int lm;
    double vk;
    double exri;
    // This vector will be passed to the new object. DO NOT DELETE HERE!
    dcomplex *elems;
    double radius = 0.0;
    status = hdf_file->read("IS", "INT32", &is);
    status = hdf_file->read("L_MAX", "INT32", &lm);
    status = hdf_file->read("VK", "FLOAT64", &vk);
    status = hdf_file->read("EXRI", "FLOAT64", &exri);
    if (is == 1111) {
      int num_elements = 2 * lm;
      double *file_vector = new double[2 * num_elements]();
      hid_t file_id = hdf_file->get_file_id();
      hid_t dset_id = H5Dopen2(file_id, "ELEMENTS", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, file_vector);
      elems = new dcomplex[num_elements]();
      for (int ei = 0; ei < num_elements; ei++) {
	elems[ei] = file_vector[2 * ei] + file_vector[2 * ei + 1] * I;
      }
      status = H5Dclose(dset_id);
      status = hdf_file->read("RADIUS", "FLOAT64", &radius);
      tm = new TransitionMatrix(is, lm, vk, exri, elems, radius);
      delete[] file_vector;
    } else if (is == 1) {
      int nlemt = 2 * lm * (lm + 2);
      int num_elements = nlemt * nlemt;
      double *file_vector = new double[2 * num_elements]();
      hid_t file_id = hdf_file->get_file_id();
      hid_t dset_id = H5Dopen2(file_id, "ELEMENTS", H5P_DEFAULT);
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, file_vector);
      elems = new dcomplex[num_elements]();
      for (int ei = 0; ei < num_elements; ei++) {
	elems[ei] = file_vector[2 * ei] + file_vector[2 * ei + 1] * I;
      }
      status = H5Dclose(dset_id);
      tm = new TransitionMatrix(is, lm, vk, exri, elems, radius);
      delete[] file_vector;
    }
    status = hdf_file->close();
  } else {
    printf("[ERROR] could not open file \"%s\"\n", file_name.c_str());
  }
  delete hdf_file;
  return tm;
}

TransitionMatrix* TransitionMatrix::from_legacy(const std::string& file_name) {
  fstream ttms;
  TransitionMatrix *tm = NULL;
  ttms.open(file_name, ios::binary | ios::in);
  if (ttms.is_open()) {
    int num_elements = 0;
    int is;
    int lm;
    double vk;
    double exri;
    // This vector will be passed to the new object. DO NOT DELETE HERE!
    dcomplex *elems;
    double radius = 0.0;
    ttms.read(reinterpret_cast<char *>(&is), sizeof(int));
    ttms.read(reinterpret_cast<char *>(&lm), sizeof(int));
    ttms.read(reinterpret_cast<char *>(&vk), sizeof(double));
    ttms.read(reinterpret_cast<char *>(&exri), sizeof(double));
    if (is == 1111) {
      num_elements = lm * 2;
      elems = new dcomplex[num_elements]();
      for (int ei = 0; ei < num_elements; ei++) {
	double vreal, vimag;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	elems[ei] = vreal + vimag * I;
      }
      ttms.read(reinterpret_cast<char *>(&radius), sizeof(double));
      tm = new TransitionMatrix(is, lm, vk, exri, elems, radius);
    } else if (is == 1) {
      int nlemt = 2 * lm * (lm + 2);
      num_elements = nlemt * nlemt;
      elems = new dcomplex[num_elements]();
      for (int ei = 0; ei < num_elements; ei++) {
	double vreal, vimag;
	ttms.read(reinterpret_cast<char *>(&vreal), sizeof(double));
	ttms.read(reinterpret_cast<char *>(&vimag), sizeof(double));
	elems[ei] = vreal + vimag * I;
      }
      tm = new TransitionMatrix(is, lm, vk, exri, elems);
    }
  } else {
    printf("[ERROR] could not open file \"%s\"\n", file_name.c_str());
  }
  return tm;
}

void TransitionMatrix::write_binary(const std::string& file_name, const std::string& mode) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"";
    throw UnrecognizedFormatException(message);
  }
}

void TransitionMatrix::write_binary(
  const std::string& file_name, np_int nlemt, int lm, double vk,
  double exri, dcomplex **am0m, const std::string& mode
) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name, nlemt, lm, vk, exri, am0m);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name, nlemt, lm, vk, exri, am0m);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"";
    throw UnrecognizedFormatException(message);
  }
}

void TransitionMatrix::write_binary(
  const std::string& file_name, int lm, double vk, double exri,
  dcomplex **rmi, dcomplex **rei, double radius, const std::string& mode
) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name, lm, vk, exri, rmi, rei, radius);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name, lm, vk, exri, rmi, rei, radius);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"";
    throw UnrecognizedFormatException(message);
  }
}

void TransitionMatrix::write_hdf5(const std::string& file_name) {
  check_output_file(file_name, "T-matrix file");
  if (_is == 1 || _is == 1111) {
    List<string> rec_name_list(1);
    List<string> rec_type_list(1);
    List<void *> rec_ptr_list(1);
    string str_type, str_name;
    rec_name_list.set(0, "IS");
    rec_type_list.set(0, "INT32_(1)");
    rec_ptr_list.set(0, &_is);
    rec_name_list.append("L_MAX");
    rec_type_list.append("INT32_(1)");
    rec_ptr_list.append(&_l_max);
    rec_name_list.append("VK");
    rec_type_list.append("FLOAT64_(1)");
    rec_ptr_list.append(&_vk);
    rec_name_list.append("EXRI");
    rec_type_list.append("FLOAT64_(1)");
    rec_ptr_list.append(&_exri);
    rec_name_list.append("ELEMENTS");
    str_type = "COMPLEX128_(" + to_string(shape[0]) + "," + to_string(shape[1]) + ")";
    rec_type_list.append(str_type);
    dcomplex *p_first = elements;
    rec_ptr_list.append(p_first);
    if (is == 1111) {
      rec_name_list.append("RADIUS");
      rec_type_list.append("FLOAT64_(1)");
      rec_ptr_list.append(&_sphere_radius);
    }

    string *rec_names = rec_name_list.to_array();
    string *rec_types = rec_type_list.to_array();
    void **rec_pointers = rec_ptr_list.to_array();
    const int rec_num = rec_name_list.length();
    FileSchema schema(rec_num, rec_types, rec_names);
    HDFFile *hdf_file = HDFFile::from_schema(schema, file_name, H5F_ACC_TRUNC);
    for (int ri = 0; ri < rec_num; ri++)
      hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
    hdf_file->close();
    
    p_first = NULL;
    delete[] rec_names;
    delete[] rec_types;
    delete[] rec_pointers;
    delete hdf_file;
  } else {
    string message = "Unrecognized matrix data.";
    throw UnrecognizedFormatException(message);
  }
}

void TransitionMatrix::write_hdf5(
  const std::string& file_name, np_int nlemt, int lm, double vk, double exri, dcomplex **am0m
) {
  check_output_file(file_name, "T-matrix file");
  int is = 1;
  List<string> rec_name_list(1);
  List<string> rec_type_list(1);
  List<void *> rec_ptr_list(1);
  string str_type, str_name;
  rec_name_list.set(0, "IS");
  rec_type_list.set(0, "INT32_(1)");
  rec_ptr_list.set(0, &is);
  rec_name_list.append("L_MAX");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&lm);
  rec_name_list.append("VK");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&vk);
  rec_name_list.append("EXRI");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&exri);
  rec_name_list.append("ELEMENTS");
  str_type = "COMPLEX128_(" + to_string(nlemt) + "," + to_string(nlemt) + ")";
  rec_type_list.append(str_type);
  // The (N x M) matrix of complex is converted to a (N x 2M) matrix of double,
  // where REAL(E_i,j) -> E_i,(2 j) and IMAG(E_i,j) -> E_i,(2 j + 1)
  dcomplex *p_first = am0m[0];
  rec_ptr_list.append(p_first);
  
  string *rec_names = rec_name_list.to_array();
  string *rec_types = rec_type_list.to_array();
  void **rec_pointers = rec_ptr_list.to_array();
  const int rec_num = rec_name_list.length();
  FileSchema schema(rec_num, rec_types, rec_names);
  HDFFile *hdf_file = HDFFile::from_schema(schema, file_name, H5F_ACC_TRUNC);
  for (int ri = 0; ri < rec_num; ri++)
    hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
  hdf_file->close();

  p_first = NULL;
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete hdf_file;
}

void TransitionMatrix::write_hdf5(
  const std::string& file_name, int lm, double vk, double exri,
  dcomplex **rmi, dcomplex **rei, double radius
) {
  check_output_file(file_name, "T-matrix file");
  int is = 1111;
  List<string> rec_name_list(1);
  List<string> rec_type_list(1);
  List<void *> rec_ptr_list(1);
  string str_type, str_name;
  rec_name_list.set(0, "IS");
  rec_type_list.set(0, "INT32_(1)");
  rec_ptr_list.set(0, &is);
  rec_name_list.append("L_MAX");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&lm);
  rec_name_list.append("VK");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&vk);
  rec_name_list.append("EXRI");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&exri);
  dcomplex *elems = new dcomplex[2 * lm]();
  for (int ei = 0; ei < lm; ei++) {
    elems[2 * ei] = -1.0 / rmi[ei][0];
    elems[2 * ei + 1] = -1.0 / rei[ei][0];
  }
  rec_name_list.append("ELEMENTS");
  str_type = "COMPLEX128_(" + to_string(lm) + "," + to_string(2) + ")";
  rec_type_list.append(str_type);
  rec_ptr_list.append(elems);
  rec_name_list.append("RADIUS");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&radius);

  string *rec_names = rec_name_list.to_array();
  string *rec_types = rec_type_list.to_array();
  void **rec_pointers = rec_ptr_list.to_array();
  const int rec_num = rec_name_list.length();
  FileSchema schema(rec_num, rec_types, rec_names);
  HDFFile *hdf_file = HDFFile::from_schema(schema, file_name, H5F_ACC_TRUNC);
  for (int ri = 0; ri < rec_num; ri++)
    hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
  hdf_file->close();
    
  delete[] elems;
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete hdf_file;
}

void TransitionMatrix::write_legacy(const std::string& file_name) {
  check_output_file(file_name, "T-matrix file");
  fstream ttms;
  if (_is == 1111 || _is == 1) {
    ttms.open(file_name, ios::binary | ios::out);
    if (ttms.is_open()) {
      ttms.write(reinterpret_cast<char *>(&_is), sizeof(int));
      ttms.write(reinterpret_cast<char *>(&_l_max), sizeof(int));
      ttms.write(reinterpret_cast<char *>(&_vk), sizeof(double));
      ttms.write(reinterpret_cast<char *>(&_exri), sizeof(double));
    }
  } else {
    string message = "Unrecognized matrix data.";
    throw UnrecognizedFormatException(message);
  }
  if (ttms.is_open()) {
    int num_elements = shape[0] * shape[1];
    for (int ei = 0; ei < num_elements; ei++) {
      dcomplex element1 = elements[ei];
      double vreal, vimag;
      vreal = real(element1);
      vimag = imag(element1);
      ttms.write(reinterpret_cast<char *>(&vreal), sizeof(double));
      ttms.write(reinterpret_cast<char *>(&vimag), sizeof(double));
    }
    if (_is == 1111) {
      ttms.write(reinterpret_cast<char *>(&_sphere_radius), sizeof(double));
    }
    ttms.close();
  } else { // Failed to open output file. Should never happen.
    throw OutputError("cannot write the T-matrix file \"" + file_name + "\"");
  }
}

void TransitionMatrix::write_legacy(
  const std::string& file_name, np_int nlemt, int lm, double vk,
  double exri, dcomplex **am0m
) {
  check_output_file(file_name, "T-matrix file");
  fstream ttms;
  int is = 1;
  ttms.open(file_name, ios::binary | ios::out);
  if (ttms.is_open()) {
    ttms.write(reinterpret_cast<char *>(&is), sizeof(int));
    ttms.write(reinterpret_cast<char *>(&lm), sizeof(int));
    ttms.write(reinterpret_cast<char *>(&vk), sizeof(double));
    ttms.write(reinterpret_cast<char *>(&exri), sizeof(double));
    double rval, ival;
    for (np_int ei = 0; ei < nlemt; ei++) {
      for (np_int ej = 0; ej < nlemt; ej++) {
	rval = real(am0m[ei][ej]);
	ival = imag(am0m[ei][ej]);
	ttms.write(reinterpret_cast<char *>(&rval), sizeof(double));
	ttms.write(reinterpret_cast<char *>(&ival), sizeof(double));
      }
    }
    ttms.close();
  } else { // Failed to open output file. Should never happen.
    throw OutputError("cannot write the T-matrix file \"" + file_name + "\"");
  }
}

void TransitionMatrix::write_legacy(
  const std::string& file_name, int lm, double vk, double exri,
  dcomplex **rmi, dcomplex **rei, double radius
) {
  check_output_file(file_name, "T-matrix file");
  fstream ttms;
  int is = 1111;
  ttms.open(file_name, ios::binary | ios::out);
  if (ttms.is_open()) {
    ttms.write(reinterpret_cast<char *>(&is), sizeof(int));
    ttms.write(reinterpret_cast<char *>(&lm), sizeof(int));
    ttms.write(reinterpret_cast<char *>(&vk), sizeof(double));
    ttms.write(reinterpret_cast<char *>(&exri), sizeof(double));
    double rval, ival;
    dcomplex element;
    for (int ei = 0; ei < lm; ei++) {
      element = -1.0 / rmi[ei][0];
      rval = real(element);
      ival = imag(element);
      ttms.write(reinterpret_cast<char *>(&rval), sizeof(double));
      ttms.write(reinterpret_cast<char *>(&ival), sizeof(double));
      element = -1.0 / rei[ei][0];
      rval = real(element);
      ival = imag(element);
      ttms.write(reinterpret_cast<char *>(&rval), sizeof(double));
      ttms.write(reinterpret_cast<char *>(&ival), sizeof(double));
    }
    ttms.write(reinterpret_cast<char *>(&radius), sizeof(double));
    ttms.close();
  } else { // Failed to open output file. Should never happen.
    throw OutputError("cannot write the T-matrix file \"" + file_name + "\"");
  }
}

bool TransitionMatrix::operator==(const TransitionMatrix &other) const {
  if (_is != other._is) {
    return false;
  }
  if (_l_max != other._l_max) {
    return false;
  }
  if (_vk != other._vk) {
    return false;
  }
  if (_exri != other._exri) {
    return false;
  }
  if (_sphere_radius != other._sphere_radius) {
    return false;
  }
  if (shape[0] != other.shape[0]) {
    return false;
  }
  if (shape[1] != other.shape[1]) {
    return false;
  }
  int num_elements = shape[0] * shape[1];
  for (int ei = 0; ei < num_elements; ei++) {
    if (elements[ei] != other.elements[ei]) {
      return false;
    }
  }
  return true;
}

bool TransitionMatrix::operator!=(const TransitionMatrix &other) const {
  return !(*this == other);
}

void TransitionMatrix::mpisend(const mixMPI *mpidata, int destination) {
#ifdef MPI_VERSION
  int ints[4] = {_is, _l_max, shape[0], shape[1]};
  double doubles[3] = {_vk, _exri, _sphere_radius};
  MPI_Send(ints, 4, MPI_INT, destination, 30, MPI_COMM_WORLD);
  MPI_Send(doubles, 3, MPI_DOUBLE, destination, 31, MPI_COMM_WORLD);
  MPI_Send(elements, shape[0] * shape[1], MPI_C_DOUBLE_COMPLEX, destination, 32, MPI_COMM_WORLD);
#endif
}

TransitionMatrix* TransitionMatrix::mpireceive(const mixMPI *mpidata, int source) {
  TransitionMatrix *tm = NULL;
#ifdef MPI_VERSION
  int ints[4];
  double doubles[3];
  MPI_Recv(ints, 4, MPI_INT, source, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(doubles, 3, MPI_DOUBLE, source, 31, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  const int num_elements = ints[2] * ints[3];
  dcomplex *elems = new dcomplex[num_elements];
  MPI_Recv(elems, num_elements, MPI_C_DOUBLE_COMPLEX, source, 32, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  tm = new TransitionMatrix(ints[0], ints[1], doubles[0], doubles[1], elems, doubles[2]);
  // The shape of the sender, as the elements were sent.
  tm->shape[0] = ints[2];
  tm->shape[1] = ints[3];
#endif
  return tm;
}

void TransitionMatrix::write_from_main_process(
  const mixMPI *mpidata, TransitionMatrix *tm, const std::string& output_path
) {
  TransitionMatrix *received = NULL;
#ifdef MPI_VERSION
  if (mpidata->mpirunning) {
    // The process that holds the T-matrix (the highest rank, if several did).
    int holder = (tm != NULL) ? mpidata->rank : -1;
    int global_holder = -1;
    MPI_Allreduce(&holder, &global_holder, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_holder > 0) {
      if (mpidata->rank == global_holder) {
	tm->mpisend(mpidata, 0);
      } else if (mpidata->rank == 0) {
	received = TransitionMatrix::mpireceive(mpidata, global_holder);
	tm = received;
      }
    }
  }
#endif
  if (mpidata->rank == 0 && tm != NULL) {
    try {
      tm->write_binary(output_path + "/c_TTMS.hd5", "HDF5");
      tm->write_binary(output_path + "/c_TTMS", "LEGACY");
    } catch (...) {
      delete received;
      throw;
    }
  }
  delete received;
}
