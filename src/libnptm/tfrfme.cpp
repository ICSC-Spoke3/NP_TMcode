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

/*! \file tfrfme.cpp
 *
 * \brief Implementation of the trapping calculation objects.
 */

#include <exception>
#include <fstream>
#include <hdf5.h>
#include <string>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_TFRFME_H_
#include "../include/tfrfme.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifdef USE_TARGET_OFFLOAD
#include <cstdlib>
#endif

using namespace std;

// >>> START OF Swap1 CLASS IMPLEMENTATION <<<
Swap1::Swap1(int lm, int nkv) {
  _nkv = nkv;
  _nlmmt = 2 * lm * (lm + 2);
  const int size = nkv * nkv * _nlmmt;
  _wk = new dcomplex[size]();
  wk = _wk;
  _last_index = 0;
}

Swap1* Swap1::from_binary(const std::string& file_name, const std::string& mode) {
  Swap1 *instance = NULL;
  if (mode.compare("LEGACY") == 0) {
    instance = from_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    instance = from_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"\n";
    throw UnrecognizedFormatException(message);
  }
  return instance;
}

Swap1* Swap1::from_hdf5(const std::string& file_name) {
  Swap1 *instance = NULL;
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(file_name, flags);
  herr_t status = hdf_file->get_status();
  double *elements;
  string str_type;
  int nlmmt, nkv, lm, num_elements, index;
  dcomplex value;
  if (status == 0) {
    status = hdf_file->read("NLMMT", "INT32", &nlmmt);
    status = hdf_file->read("NKV", "INT32", &nkv);
    lm = (int)(sqrt(4.0 + 2.0 * nlmmt) / 2.0) - 1;
    num_elements = 2 * nlmmt * nkv * nkv;
    instance = new Swap1(lm, nkv);
    elements = new double[num_elements]();
    str_type = "FLOAT64_(" + to_string(num_elements) + ")";
    status = hdf_file->read("WK", str_type, elements);
    for (int wi = 0; wi < num_elements / 2; wi++) {
      index = 2 * wi;
      value = elements[index] + elements[index + 1] * I;
      instance->_wk[wi] = value;
    } // wi loop
    delete[] elements;
    status = hdf_file->close();
    delete hdf_file;
  }
  return instance;
}

Swap1* Swap1::from_legacy(const std::string& file_name) {
  fstream input;
  Swap1 *instance = NULL;
  int nlmmt, nkv, lm;
  double rval, ival;
  input.open(file_name.c_str(), ios::in | ios::binary);
  if (input.is_open()) {
    input.read(reinterpret_cast<char *>(&nlmmt), sizeof(int));
    lm = (int)(sqrt(4.0 + 2.0 * nlmmt) / 2.0) - 1;
    input.read(reinterpret_cast<char *>(&nkv), sizeof(int));
    instance = new Swap1(lm, nkv);
    int num_elements = nlmmt * nkv * nkv;
    for (int j = 0; j < num_elements; j++) {
      input.read(reinterpret_cast<char *>(&rval), sizeof(double));
      input.read(reinterpret_cast<char *>(&ival), sizeof(double));
      instance->_wk[j] = rval + ival * I;
    }
    input.close();
  } else {
    printf("ERROR: could not open input file \"%s\"\n", file_name.c_str());
  }
  return instance;
}

long Swap1::get_size(int lm, int nkv) {
  long size = (long)(3 * sizeof(int));
  size += (long)(sizeof(dcomplex) * 2 * lm * (lm + 2) * nkv * nkv);
  return size;
}

void Swap1::write_binary(const std::string& file_name, const std::string& mode) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"\n";
    throw UnrecognizedFormatException(message);
  }
}

void Swap1::write_hdf5(const std::string& file_name) {
  List<string> rec_name_list(1);
  List<string> rec_type_list(1);
  List<void *> rec_ptr_list(1);
  herr_t status;
  string str_type;
  int num_elements = 2 * _nlmmt * _nkv * _nkv;
  rec_name_list.set(0, "NLMMT");
  rec_type_list.set(0, "INT32_(1)");
  rec_ptr_list.set(0, &_nlmmt);
  rec_name_list.append("NKV");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nkv);
  rec_name_list.append("WK");
  str_type = "FLOAT64_(" + to_string(num_elements) + ")";
  rec_type_list.append(str_type);
  double *ptr_elements = new double[num_elements]();
  for (int wi = 0; wi < num_elements / 2; wi++) {
    ptr_elements[2 * wi] = real(_wk[wi]);
    ptr_elements[2 * wi + 1] = imag(_wk[wi]);
  }
  rec_ptr_list.append(ptr_elements);

  string *rec_names = rec_name_list.to_array();
  string *rec_types = rec_type_list.to_array();
  void **rec_pointers = rec_ptr_list.to_array();
  const int rec_num = rec_name_list.length();
  FileSchema schema(rec_num, rec_types, rec_names);
  HDFFile *hdf_file = HDFFile::from_schema(schema, file_name, H5F_ACC_TRUNC);
  for (int ri = 0; ri < rec_num; ri++)
    status = hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
  status = hdf_file->close();

  delete[] ptr_elements;
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete hdf_file;
}

void Swap1::write_legacy(const std::string& file_name) {
  fstream output;
  double rval, ival;
  output.open(file_name.c_str(), ios::out | ios::binary);
  if (output.is_open()) {
    int num_elements = _nlmmt * _nkv * _nkv;
    output.write(reinterpret_cast<char *>(&_nlmmt), sizeof(int));
    output.write(reinterpret_cast<char *>(&_nkv), sizeof(int));
    for (int j = 0; j < num_elements; j++) {
      rval = real(_wk[j]);
      ival = imag(_wk[j]);
      output.write(reinterpret_cast<char *>(&rval), sizeof(double));
      output.write(reinterpret_cast<char *>(&ival), sizeof(double));
    }
    output.close();
  } else { // Should never occur
    printf("ERROR: could not open output file \"%s\"\n", file_name.c_str());
  }
}

bool Swap1::operator ==(Swap1 &other) {
  if (_nlmmt != other._nlmmt) {
    return false;
  }
  if (_nkv != other._nkv) {
    return false;
  }
  const int num_elements = _nlmmt * _nkv * _nkv;
  for (int i = 0; i < num_elements; i++) {
    if (_wk[i] != other._wk[i]) {
      return false;
    }
  }
  return true;
}
// >>> END OF Swap1 CLASS IMPLEMENTATION <<<

// >>> START OF Swap2 CLASS IMPLEMENTATION <<<
Swap2::Swap2(int nkv) {
  _nkv = nkv;
#ifdef USE_TARGET_OFFLOAD
  vkv = (double *)aligned_alloc(64, _nkv * sizeof(double));
  vec_vkzm = (double *)aligned_alloc(64, _nkv * _nkv * sizeof(double));
#pragma omp parallel for collapse(2)
  for (int i = 0; i < _nkv; i++) {
    for (int j = 0; j < _nkv; j++) {
      vkv[i] = 0.0;
      vec_vkzm[_nkv * i +j] = 0.0;
    }
  }
#else
  vkv = new double[_nkv]();
  vec_vkzm = new double[_nkv * _nkv]();
#endif // USE TARGET_OFFLOAD
  _last_vector = 0;
  _last_matrix = 0;
}

Swap2::~Swap2() {
#ifdef USE_TARGET_OFFLOAD
  free(vkv);
  free(vec_vkzm);
#else
  delete[] vkv;
  delete[] vec_vkzm;
#endif // USE_TARGET_OFFLOAD
}

Swap2* Swap2::from_binary(const std::string& file_name, const std::string& mode) {
  Swap2 *instance = NULL;
  if (mode.compare("LEGACY") == 0) {
    instance = from_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    instance = from_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"\n";
    throw UnrecognizedFormatException(message);
  }
  return instance;
}

Swap2* Swap2::from_hdf5(const std::string& file_name) {
  Swap2 *instance = NULL;
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(file_name, flags);
  herr_t status = hdf_file->get_status();
  string str_type;
  int fnkv, fnlmmt, fnrvc;
  double value;
  if (status == 0) {
    status = hdf_file->read("NKV", "INT32", &fnkv);
    instance = new Swap2(fnkv);
    str_type = "FLOAT64_(" + to_string(fnkv) + ")";
    status = hdf_file->read("VKV", str_type, instance->vkv);
    str_type = "FLOAT64_(" + to_string(fnkv * fnkv) + ")";
    status = hdf_file->read("VKZM", str_type, instance->vec_vkzm);
    status = hdf_file->read("APFAFA", "FLOAT64", &value);
    instance->set_param("apfafa", value);
    status = hdf_file->read("PMF", "FLOAT64", &value);
    instance->set_param("pmf", value);
    status = hdf_file->read("SPD", "FLOAT64", &value);
    instance->set_param("spd", value);
    status = hdf_file->read("RIR", "FLOAT64", &value);
    instance->set_param("rir", value);
    status = hdf_file->read("FTCN", "FLOAT64", &value);
    instance->set_param("ftcn", value);
    status = hdf_file->read("FSHMX", "FLOAT64", &value);
    instance->set_param("fshmx", value);
    status = hdf_file->read("VXYZMX", "FLOAT64", &value);
    instance->set_param("vxyzmx", value);
    status = hdf_file->read("DELXYZ", "FLOAT64", &value);
    instance->set_param("delxyz", value);
    status = hdf_file->read("VKNMX", "FLOAT64", &value);
    instance->set_param("vknmx", value);
    status = hdf_file->read("DELK", "FLOAT64", &value);
    instance->set_param("delk", value);
    status = hdf_file->read("DELKS", "FLOAT64", &value);
    instance->set_param("delks", value);
    status = hdf_file->read("NLMMT", "INT32", &fnlmmt);
    instance->set_param("nlmmt", 1.0 * fnlmmt);
    status = hdf_file->read("NRVC", "INT32", &fnrvc);
    instance->set_param("nrvc", 1.0 * fnrvc);
    status = hdf_file->close();
    delete hdf_file;
  }
  return instance;
}

Swap2* Swap2::from_legacy(const std::string& file_name) {
  fstream input;
  Swap2 *instance = NULL;
  int fnkv, fnlmmt, fnrvc;
  double *fvkzm = NULL;
  double *fvkv = NULL;
  double value;
  input.open(file_name.c_str(), ios::in | ios::binary);
  if (input.is_open()) {
    input.read(reinterpret_cast<char *>(&fnkv), sizeof(int));
    instance = new Swap2(fnkv);
    fvkzm = instance->vec_vkzm;
    fvkv = instance->get_vector();
    for (int vj = 0; vj < fnkv; vj++) {
      input.read(reinterpret_cast<char *>(&value), sizeof(double));
      fvkv[vj] = value;
    }
    for (int mi = 0; mi < fnkv; mi++) {
      for (int mj = 0; mj < fnkv; mj++) {
	input.read(reinterpret_cast<char *>(&value), sizeof(double));
	fvkzm[fnkv * mi + mj] = value;
      }
    }
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("apfafa", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("pmf", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("spd", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("rir", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("ftcn", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("fshmx", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("vxyzmx", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("delxyz", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("vknmx", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("delk", value);
    input.read(reinterpret_cast<char *>(&value), sizeof(double));
    instance->set_param("delks", value);
    input.read(reinterpret_cast<char *>(&fnlmmt), sizeof(int));
    instance->set_param("nlmmt", 1.0 * fnlmmt);
    input.read(reinterpret_cast<char *>(&fnrvc), sizeof(int));
    instance->set_param("nrvc", 1.0 * fnrvc);
    input.close();
  } else {
    printf("ERROR: could not open input file \"%s\"\n", file_name.c_str());
  }
  return instance;
}

long Swap2::get_size(int nkv) {
  long size = (long)(10 * sizeof(int) + 22 * sizeof(double));
  size += (long)(sizeof(double) * nkv * (nkv + 1));
  return size;
}

void Swap2::push_matrix(double value) {
  int col = _last_matrix % (_nkv - 1);
  int row = _last_matrix - (_nkv * row);
  vec_vkzm[nkv * row + col] = value;
  _last_matrix++;
}

void Swap2::set_param(const std::string& param_name, double value) {
  if (param_name.compare("nkv") == 0) _nkv = (int)value;
  else if (param_name.compare("apfafa") == 0) _apfafa = value;
  else if (param_name.compare("pmf") == 0) _pmf = value;
  else if (param_name.compare("spd") == 0) _spd = value;
  else if (param_name.compare("rir") == 0) _rir = value;
  else if (param_name.compare("ftcn") == 0) _ftcn = value;
  else if (param_name.compare("fshmx") == 0) _fshmx = value;
  else if (param_name.compare("vxyzmx") == 0) _vxyzmx = value;
  else if (param_name.compare("delxyz") == 0) _delxyz = value;
  else if (param_name.compare("vknmx") == 0) _vknmx = value;
  else if (param_name.compare("delk") == 0) _delk = value;
  else if (param_name.compare("delks") == 0) _delks = value;
  else if (param_name.compare("nlmmt") == 0) _nlmmt = (int)value;
  else if (param_name.compare("nrvc") == 0) _nrvc = (int)value;
  else {
    string message = "Unrecognized parameter name \"" + param_name + "\"\n";
    throw UnrecognizedParameterException(message);
  }
}

void Swap2::write_binary(const std::string& file_name, const std::string& mode) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"\n";
    throw UnrecognizedFormatException(message);
  }
}

void Swap2::write_hdf5(const std::string& file_name) {
  List<string> rec_name_list(1);
  List<string> rec_type_list(1);
  List<void *> rec_ptr_list(1);
  herr_t status;
  string str_type;
  rec_name_list.set(0, "NKV");
  rec_type_list.set(0, "INT32_(1)");
  rec_ptr_list.set(0, &_nkv);
  rec_name_list.append("VKV");
  str_type = "FLOAT64_(" + to_string(_nkv) + ")";
  rec_type_list.append(str_type);
  rec_ptr_list.append(vkv);
  rec_name_list.append("VKZM");
  str_type = "FLOAT64_(" + to_string(_nkv * _nkv) + ")";
  rec_type_list.append(str_type);
  rec_ptr_list.append(vec_vkzm);
  rec_name_list.append("APFAFA");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_apfafa);
  rec_name_list.append("PMF");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_pmf);
  rec_name_list.append("SPD");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_spd);
  rec_name_list.append("RIR");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_rir);
  rec_name_list.append("FTCN");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_ftcn);
  rec_name_list.append("FSHMX");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_fshmx);
  rec_name_list.append("VXYZMX");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_vxyzmx);
  rec_name_list.append("delxyz");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_delxyz);
  rec_name_list.append("VKNMX");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_vknmx);
  rec_name_list.append("DELK");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_delk);
  rec_name_list.append("DELKS");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_delks);
  rec_name_list.append("NLMMT");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nlmmt);
  rec_name_list.append("NRVC");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nrvc);

  string *rec_names = rec_name_list.to_array();
  string *rec_types = rec_type_list.to_array();
  void **rec_pointers = rec_ptr_list.to_array();
  const int rec_num = rec_name_list.length();
  FileSchema schema(rec_num, rec_types, rec_names);
  HDFFile *hdf_file = HDFFile::from_schema(schema, file_name, H5F_ACC_TRUNC);
  for (int ri = 0; ri < rec_num; ri++)
    status = hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
  status = hdf_file->close();

  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete hdf_file;
}

void Swap2::write_legacy(const std::string& file_name) {
  fstream output;
  double value;
  output.open(file_name.c_str(), ios::out | ios::binary);
  if (output.is_open()) {
    output.write(reinterpret_cast<char *>(&_nkv), sizeof(int));
    for (int j = 0; j < _nkv; j++){
      value = vkv[j];
      output.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }
    for (int mi = 0; mi < _nkv; mi++) {
      for (int mj = 0; mj < _nkv; mj++) {
	value = vec_vkzm[nkv * mi + mj];
	output.write(reinterpret_cast<const char*>(&value), sizeof(double));
      }
    }
    output.write(reinterpret_cast<const char*>(&_apfafa), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_pmf), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_spd), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_rir), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_ftcn), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_fshmx), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_vxyzmx), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_delxyz), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_vknmx), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_delk), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_delks), sizeof(double));
    output.write(reinterpret_cast<const char*>(&_nlmmt), sizeof(int));
    output.write(reinterpret_cast<const char*>(&_nrvc), sizeof(int));
    output.close();
  } else { // Should never occur
    printf("ERROR: could not open output file \"%s\"\n", file_name.c_str());
  }
}

bool Swap2::operator ==(Swap2 &other) {
  if (_nlmmt != other._nlmmt) {
    return false;
  }
  if (_nrvc != other._nrvc) {
    return false;
  }
  if (_nkv != other._nkv) {
    return false;
  }
  if (_apfafa != other._apfafa) {
    return false;
  }
  if (_pmf != other._pmf) {
    return false;
  }
  if (_spd != other._spd) {
    return false;
  }
  if (_rir != other._rir) {
    return false;
  }
  if (_ftcn != other._ftcn) {
    return false;
  }
  if (_fshmx != other._fshmx) {
    return false;
  }
  if (_vxyzmx != other._vxyzmx) {
    return false;
  }
  if (_delxyz != other._delxyz) {
    return false;
  }
  if (_vknmx != other._vknmx) {
    return false;
  }
  if (_delk != other._delk) {
    return false;
  }
  if (_delks != other._delks) {
    return false;
  }
  for (int vi = 0; vi < _nkv; vi++) {
    if (vkv[vi] != other.vkv[vi]) {
      return false;
    }
  }
  for (int mi = 0; mi < _nkv; mi++) {
    int nkvi = nkv * mi;
    for (int mj = 0; mj < _nkv; mj++) {
      if (vec_vkzm[nkvi + mj] != other.vec_vkzm[nkvi + mj]) {
	return false;
      }
    }
  }
  return true;
}
// >>> END OF Swap2 CLASS IMPLEMENTATION <<<

// >>> START OF TFRFME CLASS IMPLEMENTATION <<<
TFRFME::TFRFME(int lmode, int lm, int nkv, int nxv, int nyv, int nzv) {
  _lmode = lmode;
  _lm = lm;
  _nkv = nkv;
  _nxv = nxv;
  _nyv = nyv;
  _nzv = nzv;
  _vk = 0.0;
  _exri = 0.0;
  _an = 0.0;
  _ff = 0.0;
  _tra = 0.0;
  _spd = 0.0;
  _frsh = 0.0;
  _exril = 0.0;

  // Array initialization
  _nlmmt = _lm * (_lm + 2) * 2;
  _nrvc = _nxv * _nyv * _nzv;
#ifdef USE_TARGET_OFFLOAD
  xv = (double *)aligned_alloc(64, sizeof(double) * _nxv);
  yv = (double *)aligned_alloc(64, sizeof(double) * _nyv);
  zv = (double *)aligned_alloc(64, sizeof(double) * _nzv);
  vec_wsum = (dcomplex *)aligned_alloc(64, sizeof(dcomplex) * _nrvc * _nlmmt);
#else
  xv = new double[_nxv]();
  yv = new double[_nyv]();
  zv = new double[_nzv]();
  vec_wsum = new dcomplex[_nrvc * _nlmmt]();
#endif // USE_TARGET_OFFLOAD
}

TFRFME::~TFRFME() {
#ifdef USE_TARGET_OFFLOAD
  free(xv);
  free(yv);
  free(zv);
  free(vec_wsum);
#else
  delete[] xv;
  delete[] yv;
  delete[] zv;
  delete[] vec_wsum;
#endif
}

TFRFME* TFRFME::from_binary(const std::string& file_name, const std::string& mode) {
  TFRFME *instance = NULL;
  if (mode.compare("LEGACY") == 0) {
    instance = from_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    instance = from_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"\n";
    throw UnrecognizedFormatException(message);
  }
  return instance;
}

TFRFME* TFRFME::from_hdf5(const std::string& file_name) {
  TFRFME *instance = NULL;
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(file_name, flags);
  herr_t status = hdf_file->get_status();
  double *elements;
  string str_type;
  int nlmmt, nrvc, num_elements;
  dcomplex value;
  if (status == 0) {
    int lmode, lm, nkv, nxv, nyv, nzv;
    double vk, exri, an, ff, tra, spd, frsh, exril;
    status = hdf_file->read("LMODE", "INT32", &lmode);
    status = hdf_file->read("LM", "INT32", &lm);
    status = hdf_file->read("NKV", "INT32", &nkv);
    status = hdf_file->read("NXV", "INT32", &nxv);
    status = hdf_file->read("NYV", "INT32", &nyv);
    status = hdf_file->read("NZV", "INT32", &nzv);
    status = hdf_file->read("VK", "FLOAT64", &vk);
    status = hdf_file->read("EXRI", "FLOAT64", &exri);
    status = hdf_file->read("AN", "FLOAT64", &an);
    status = hdf_file->read("FF", "FLOAT64", &ff);
    status = hdf_file->read("TRA", "FLOAT64", &tra);
    status = hdf_file->read("SPD", "FLOAT64", &spd);
    status = hdf_file->read("FRSH", "FLOAT64", &frsh);
    status = hdf_file->read("EXRIL", "FLOAT64", &exril);
    instance = new TFRFME(lmode, lm, nkv, nxv, nyv, nzv);
    instance->set_param("vk", vk);
    instance->set_param("exri", exri);
    instance->set_param("an", an);
    instance->set_param("ff", ff);
    instance->set_param("tra", tra);
    instance->set_param("spd", spd);
    instance->set_param("frsh", frsh);
    instance->set_param("exril", exril);
    str_type = "FLOAT64_(" + to_string(nxv) + ")";
    status = hdf_file->read("XVEC", str_type, instance->xv);
    str_type = "FLOAT64_(" + to_string(nyv) + ")";
    status = hdf_file->read("YVEC", str_type, instance->yv);
    str_type = "FLOAT64_(" + to_string(nzv) + ")";
    status = hdf_file->read("ZVEC", str_type, instance->zv);
    nlmmt = 2 * lm * (lm + 2);
    nrvc = nxv * nyv * nzv;
    num_elements = 2 * nlmmt * nrvc;
    elements = new double[num_elements]();
    str_type = "FLOAT64_(" + to_string(num_elements) + ")";
    status = hdf_file->read("WSUM", str_type, elements);
    int index = 0;
    for (int wj = 0; wj < nrvc; wj++) {
      for (int wi = 0; wi < nlmmt; wi++) {
	value = elements[index] + elements[index + 1] * I;
	instance->vec_wsum[nrvc * wi + wj] = value;
	index += 2;
      } // wi loop
    } // wj loop
    delete[] elements;
    status = hdf_file->close();
    delete hdf_file;
  }
  return instance;
}

TFRFME* TFRFME::from_legacy(const std::string& file_name) {
  fstream input;
  TFRFME *instance = NULL;
  double *coord_vec = NULL;
  input.open(file_name.c_str(), ios::in | ios::binary);
  if (input.is_open()) {
    int lmode, lm, nkv, nxv, nyv, nzv;
    double vk, exri, an, ff, tra, spd, frsh, exril;
    double dval, rval, ival;
    input.read(reinterpret_cast<char *>(&lmode), sizeof(int));
    input.read(reinterpret_cast<char *>(&lm), sizeof(int));
    input.read(reinterpret_cast<char *>(&nkv), sizeof(int));
    input.read(reinterpret_cast<char *>(&nxv), sizeof(int));
    input.read(reinterpret_cast<char *>(&nyv), sizeof(int));
    input.read(reinterpret_cast<char *>(&nzv), sizeof(int));
    input.read(reinterpret_cast<char *>(&vk), sizeof(double));
    input.read(reinterpret_cast<char *>(&exri), sizeof(double));
    input.read(reinterpret_cast<char *>(&an), sizeof(double));
    input.read(reinterpret_cast<char *>(&ff), sizeof(double));
    input.read(reinterpret_cast<char *>(&tra), sizeof(double));
    input.read(reinterpret_cast<char *>(&spd), sizeof(double));
    input.read(reinterpret_cast<char *>(&frsh), sizeof(double));
    input.read(reinterpret_cast<char *>(&exril), sizeof(double));
    instance = new TFRFME(lmode, lm, nkv, nxv, nyv, nzv);
    instance->set_param("vk", vk);
    instance->set_param("exri", exri);
    instance->set_param("an", an);
    instance->set_param("ff", ff);
    instance->set_param("tra", tra);
    instance->set_param("spd", spd);
    instance->set_param("frsh", frsh);
    instance->set_param("exril", exril);
    coord_vec = instance->get_x();
    for (int xi = 0; xi < nxv; xi++) {
      input.read(reinterpret_cast<char *>(&dval), sizeof(double));
      coord_vec[xi] = dval;
    }
    coord_vec = instance->get_y();
    for (int yi = 0; yi < nyv; yi++) {
      input.read(reinterpret_cast<char *>(&dval), sizeof(double));
      coord_vec[yi] = dval;
    }
    coord_vec = instance->get_z();
    for (int zi = 0; zi < nzv; zi++) {
      input.read(reinterpret_cast<char *>(&dval), sizeof(double));
      coord_vec[zi] = dval;
    }
    int nlmmt = 2 * lm * (lm + 2);
    int nrvc = nxv * nyv * nzv;
    for (int wj = 0; wj < nrvc; wj++) {
      for (int wi = 0; wi < nlmmt; wi++) {
	input.read(reinterpret_cast<char *>(&rval), sizeof(double));
	input.read(reinterpret_cast<char *>(&ival), sizeof(double));
	dcomplex value = rval + ival * I;
	instance->vec_wsum[nrvc * wi + wj] = value;
      } // wi loop
    } // wj loop
    input.close();
  } else {
    printf("ERROR: could not open input file \"%s\"\n", file_name.c_str());
  }
  return instance;
}

long TFRFME::get_size(int lm, int nkv, int nxv, int nyv, int nzv) {
  int nlmmt = 2 * lm * (lm + 2);
  int nrvc = nxv * nyv * nzv;
  long size = (long)(16 * sizeof(int) + 16 * sizeof(double));
  size += (long)((nxv + nyv + nzv) * sizeof(double));
  size += (long)(nlmmt * nrvc * sizeof(dcomplex));
  return size;
}

void TFRFME::set_param(const std::string& param_name, double value) {
  if (param_name.compare("vk") == 0) _vk = value;
  else if (param_name.compare("exri") == 0) _exri = value;
  else if (param_name.compare("an") == 0) _an = value;
  else if (param_name.compare("ff") == 0) _ff = value;
  else if (param_name.compare("tra") == 0) _tra = value;
  else if (param_name.compare("spd") == 0) _spd = value;
  else if (param_name.compare("frsh") == 0) _frsh = value;
  else if (param_name.compare("exril") == 0) _exril = value;
  else {
    string message = "Unrecognized parameter name \"" + param_name + "\"\n";
    throw UnrecognizedParameterException(message);
  }
}

void TFRFME::write_binary(const std::string& file_name, const std::string& mode) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"\n";
    throw UnrecognizedFormatException(message);
  }
}

void TFRFME::write_hdf5(const std::string& file_name) {
  List<string> rec_name_list(1);
  List<string> rec_type_list(1);
  List<void *> rec_ptr_list(1);
  herr_t status;
  string str_type;
  rec_name_list.set(0, "LMODE");
  rec_type_list.set(0, "INT32_(1)");
  rec_ptr_list.set(0, &_lmode);
  rec_name_list.append("LM");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_lm);
  rec_name_list.append("NKV");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nkv);
  rec_name_list.append("NXV");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nxv);
  rec_name_list.append("NYV");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nyv);
  rec_name_list.append("NZV");
  rec_type_list.append("INT32_(1)");
  rec_ptr_list.append(&_nzv);
  rec_name_list.append("VK");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_vk);
  rec_name_list.append("EXRI");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_exri);
  rec_name_list.append("AN");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_an);
  rec_name_list.append("FF");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_ff);
  rec_name_list.append("TRA");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_tra);
  rec_name_list.append("SPD");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_spd);
  rec_name_list.append("FRSH");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_frsh);
  rec_name_list.append("EXRIL");
  rec_type_list.append("FLOAT64_(1)");
  rec_ptr_list.append(&_exril);
  str_type = "FLOAT64_(" + to_string(nxv) + ")";
  rec_name_list.append("XVEC");
  rec_type_list.append(str_type);
  rec_ptr_list.append(xv);
  str_type = "FLOAT64_(" + to_string(nyv) + ")";
  rec_name_list.append("YVEC");
  rec_type_list.append(str_type);
  rec_ptr_list.append(yv);
  str_type = "FLOAT64_(" + to_string(nzv) + ")";
  rec_name_list.append("ZVEC");
  rec_type_list.append(str_type);
  rec_ptr_list.append(zv);
  rec_name_list.append("WSUM");
  int num_elements = 2 * nlmmt * nrvc;
  str_type = "FLOAT64_(" + to_string(num_elements) + ")";
  rec_type_list.append(str_type);
  // The (N x M) matrix of complex is converted to a vector of double
  // with length (2 * N * M)
  double *ptr_elements = new double[num_elements]();
  int index = 0;
  for (int wj = 0; wj < nrvc; wj++) {
    for (int wi = 0; wi < nlmmt; wi++) {
      ptr_elements[index++] = real(vec_wsum[nrvc * wi + wj]);
      ptr_elements[index++] = imag(vec_wsum[nrvc * wi + wj]);
    } // wi loop
  } // wj loop
  rec_ptr_list.append(ptr_elements);

  string *rec_names = rec_name_list.to_array();
  string *rec_types = rec_type_list.to_array();
  void **rec_pointers = rec_ptr_list.to_array();
  const int rec_num = rec_name_list.length();
  FileSchema schema(rec_num, rec_types, rec_names);
  HDFFile *hdf_file = HDFFile::from_schema(schema, file_name, H5F_ACC_TRUNC);
  for (int ri = 0; ri < rec_num; ri++)
    status = hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
  status = hdf_file->close();

  delete[] ptr_elements;
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete hdf_file;
}

void TFRFME::write_legacy(const std::string& file_name) {
  fstream output;
  output.open(file_name.c_str(), ios::out | ios::binary);
  if (output.is_open()) {
    output.write(reinterpret_cast<char *>(&_lmode), sizeof(int));
    output.write(reinterpret_cast<char *>(&_lm), sizeof(int));
    output.write(reinterpret_cast<char *>(&_nkv), sizeof(int));
    output.write(reinterpret_cast<char *>(&_nxv), sizeof(int));
    output.write(reinterpret_cast<char *>(&_nyv), sizeof(int));
    output.write(reinterpret_cast<char *>(&_nzv), sizeof(int));
    output.write(reinterpret_cast<char *>(&_vk), sizeof(double));
    output.write(reinterpret_cast<char *>(&_exri), sizeof(double));
    output.write(reinterpret_cast<char *>(&_an), sizeof(double));
    output.write(reinterpret_cast<char *>(&_ff), sizeof(double));
    output.write(reinterpret_cast<char *>(&_tra), sizeof(double));
    output.write(reinterpret_cast<char *>(&_spd), sizeof(double));
    output.write(reinterpret_cast<char *>(&_frsh), sizeof(double));
    output.write(reinterpret_cast<char *>(&_exril), sizeof(double));
    for (int xi = 0; xi < _nxv; xi++)
      output.write(reinterpret_cast<char *>(&(xv[xi])), sizeof(double));
    for (int yi = 0; yi < _nyv; yi++)
      output.write(reinterpret_cast<char *>(&(yv[yi])), sizeof(double));
    for (int zi = 0; zi < _nzv; zi++)
      output.write(reinterpret_cast<char *>(&(zv[zi])), sizeof(double));
    for (int wj = 0; wj < _nrvc; wj++) {
      for (int wi = 0; wi < _nlmmt; wi++) {
	double rval = real(vec_wsum[nrvc * wi + wj]);
	double ival = imag(vec_wsum[nrvc * wi + wj]);
	output.write(reinterpret_cast<char *>(&rval), sizeof(double));
	output.write(reinterpret_cast<char *>(&ival), sizeof(double));
      } // wi loop
    } // wj loop
    output.close();
  } else { // Should never occur
    printf("ERROR: could not open output file \"%s\"\n", file_name.c_str());
  }
}

bool TFRFME::operator ==(const TFRFME& other) {
  if (_lmode != other._lmode) {
    return false;
  }
  if (_lm != other._lm) {
    return false;
  }
  if (_nkv != other._nkv) {
    return false;
  }
  if (_nxv != other._nxv) {
    return false;
  }
  if (_nyv != other._nyv) {
    return false;
  }
  if (_nzv != other._nzv) {
    return false;
  }
  if (_vk != other._vk) {
    return false;
  }
  if (_exri != other._exri) {
    return false;
  }
  if (_an != other._an) {
    return false;
  }
  if (_ff != other._ff) {
    return false;
  }
  if (_tra != other._tra) {
    return false;
  }
  if (_spd != other._spd) {
    return false;
  }
  if (_frsh != other._frsh) {
    return false;
  }
  if (_exril != other._exril) {
    return false;
  }
  for (int xi = 0; xi < _nxv; xi++) {
    if (xv[xi] != other.xv[xi]) {
      return false;
    }
  }
  for (int yi = 0; yi < _nyv; yi++) {
    if (yv[yi] != other.yv[yi]) {
      return false;
    }
  }
  for (int zi = 0; zi < _nzv; zi++) {
    if (zv[zi] != other.zv[zi]) {
      return false;
    }
  }
  for (int wi = 0; wi < _nlmmt; wi++) {
    int i = _nrvc * wi;
    for (int wj = 0; wj < _nrvc; wj++) {
      if (vec_wsum[i + wj] != other.vec_wsum[i + wj]) {
	return false;
      }
    } // wj loop
  } // wi loop
  return true;
}
// >>> END OF TFRFME CLASS IMPLEMENTATION <<<
