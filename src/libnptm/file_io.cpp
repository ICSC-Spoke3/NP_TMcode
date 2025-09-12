/*! \file file_io.cpp
 *
 * \brief Implementation of file I/O operations.
 */
#include <exception>
#include <fstream>
#include <regex>
#include <cstring>
#include <string>
#include <hdf5.h>

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

/* >>> FileSchema class implementation <<< */
FileSchema::FileSchema(int num_rec, const std::string *rec_types, const std::string *rec_names) {
  num_records = num_rec;
  record_types = new string[num_rec];
  record_names = new string[num_rec];
  for (int i = 0; i < num_rec; i++) {
    record_types[i] = rec_types[i];
    if (rec_names != NULL) record_names[i] = rec_names[i];
    else record_names[i] = "/dset" + to_string(i);
  }
}

FileSchema::~FileSchema() {
  delete[] record_names;
  delete[] record_types;
}

string* FileSchema::get_record_names() {
  string *rec_names = new string[num_records];
  for (int i = 0; i < num_records; i++) rec_names[i] = record_names[i];
  return rec_names;
}

string* FileSchema::get_record_types() {
  string *rec_types = new string[num_records];
  for (int i = 0; i < num_records; i++) rec_types[i] = record_types[i];
  return rec_types;
}
/* >>> End of FileSchema class implementation <<< */

/* >>> HDFFile class implementation <<< */
HDFFile::HDFFile(const std::string& name, unsigned int flags, hid_t fcpl_id, hid_t fapl_id) {
  file_name = name;
  if (flags == H5F_ACC_EXCL || flags == H5F_ACC_TRUNC)
    file_id = H5Fcreate(name.c_str(), flags, fcpl_id, fapl_id);
  else if (flags == H5F_ACC_RDONLY || flags == H5F_ACC_RDWR)
    file_id = H5Fopen(name.c_str(), flags, fapl_id);
  id_list = new List<hid_t>(1);
  id_list->set(0, file_id);
  if (file_id != H5I_INVALID_HID) file_open_flag = true;
  status = (herr_t)0;
}

HDFFile::~HDFFile() {
  if (H5Iis_valid(file_id) > 0) status = H5Fclose(file_id);
  delete id_list;
}

herr_t HDFFile::close() {
  status = H5Fclose(file_id);
  if (status == 0) file_open_flag = false;
  return status;
}

HDFFile* HDFFile::from_schema(
			      FileSchema &schema, const std::string& name, unsigned int flags,
			      hid_t fcpl_id, hid_t fapl_id
) {
  HDFFile *hdf_file = new HDFFile(name, flags, fcpl_id, fapl_id);
  hid_t file_id = hdf_file->get_file_id();
  herr_t status;
  string *rec_types = schema.get_record_types();
  string *rec_names = schema.get_record_names();
  string known_types[] = {"INT32", "FLOAT64", "COMPLEX128", "INT64", "INT16"};
  int rec_num = schema.get_record_number();
  regex re;
  smatch m;
  for (int ri = 0; ri < rec_num; ri++) {
    int rank = 0;
    hsize_t *dims, *max_dims;
    hid_t data_type;
    string str_target = rec_types[ri];
    int type_index = 0;
    bool found_type = false;
    while (!found_type) {
      re = regex(known_types[type_index++]);
      if (regex_search(str_target, m, re)) {
	found_type = true;
	str_target = m.suffix().str();
	if (type_index == 1) data_type = H5Tcopy(H5T_NATIVE_INT);
	else if (type_index == 2) data_type = H5Tcopy(H5T_NATIVE_DOUBLE);
	else if (type_index == 3) data_type = H5Tcopy(H5T_NATIVE_DOUBLE);
	else if (type_index == 4) data_type = H5Tcopy(H5T_NATIVE_LONG);
	else if (type_index == 5) data_type = H5Tcopy(H5T_NATIVE_SHORT);
      }
      if (type_index == 5) break;
    }
    if (found_type) {
      re = regex("[0-9]+");
      string old_target = str_target;
      while (regex_search(str_target, m, re)) {
	rank++;
	str_target = m.suffix().str();
      }
      dims = new hsize_t[rank]();
      max_dims = new hsize_t[rank]();
      str_target = old_target;
      for (int ti = 0; ti < rank; ti++) {
	regex_search(str_target, m, re);
	hsize_t dim = (hsize_t)stoi(m.str());
	dims[ti] = dim;
	max_dims[ti] = dim;
	str_target = m.suffix().str();
      }
      int multiplier = (type_index == 3) ? 2 : 1;
      dims[rank - 1] *= multiplier;
      max_dims[rank - 1] *= multiplier;
      hid_t dataspace_id = H5Screate_simple(rank, dims, max_dims);
      hid_t dataset_id = H5Dcreate(
				   file_id, rec_names[ri].c_str(), data_type, dataspace_id, H5P_DEFAULT,
				   H5P_DEFAULT, H5P_DEFAULT
				   );
      status = H5Sclose(dataspace_id);
      status = H5Dclose(dataset_id);
      delete[] dims;
      delete[] max_dims;
    } else {
      string message = "unrecognized type \"" + rec_types[ri] + "\"\n";
      throw UnrecognizedParameterException(message);
    }
  }

  delete[] rec_types;
  delete[] rec_names;
  return hdf_file;
}

herr_t HDFFile::read(
		     const std::string& dataset_name, const std::string& data_type, void *buffer,
		     hid_t mem_space_id, hid_t file_space_id, hid_t dapl_id,
		     hid_t dxpl_id
) {
  string known_types[] = {"INT32", "FLOAT64", "COMPLEX128", "INT64", "INT16"};
  regex re;
  smatch m;
  bool found_type = false;
  int type_index = 0;
  while (!found_type) {
    re = regex(known_types[type_index++]);
    found_type = regex_search(data_type, m, re);
    if (type_index == 5) break;
  }
  if (found_type) {
    hid_t dataset_id = H5Dopen2(file_id, dataset_name.c_str(), dapl_id);
    hid_t mem_type_id;
    switch (type_index) {
    case 1:
      mem_type_id = H5T_NATIVE_INT; break;
    case 2:
      mem_type_id = H5T_NATIVE_DOUBLE; break;
    case 3:
      mem_type_id = H5T_NATIVE_DOUBLE; break;
    case 4:
      mem_type_id = H5T_NATIVE_LONG; break;
    case 5:
      mem_type_id = H5T_NATIVE_SHORT; break;
    default:
      throw UnrecognizedParameterException("unrecognized data type \"" + data_type + "\"");
    }
    if (dataset_id != H5I_INVALID_HID) {
      status = H5Dread(dataset_id, mem_type_id, mem_space_id, file_space_id, dxpl_id, buffer);
      if (status == 0) status = H5Dclose(dataset_id);
      else status = (herr_t)-2;
    } else {
      status = (herr_t)-1;
    }
  } else {
    throw UnrecognizedParameterException("unrecognized data type \"" + data_type + "\"");
  }
  return status;
}

herr_t HDFFile::write(
		      const std::string& dataset_name, const std::string& data_type, const void *buffer,
		      hid_t mem_space_id, hid_t file_space_id, hid_t dapl_id,
		      hid_t dxpl_id
) {
  string known_types[] = {"INT32", "FLOAT64", "COMPLEX128", "INT64", "INT16"};
  regex re;
  smatch m;
  bool found_type = false;
  int type_index = 0;
  while (!found_type) {
    re = regex(known_types[type_index++]);
    found_type = regex_search(data_type, m, re);
    if (type_index == 5) break;
  }
  if (found_type) {
    hid_t dataset_id = H5Dopen2(file_id, dataset_name.c_str(), dapl_id);
    hid_t mem_type_id;
    switch (type_index) {
    case 1:
      mem_type_id = H5T_NATIVE_INT; break;
    case 2:
      mem_type_id = H5T_NATIVE_DOUBLE; break;
    case 3:
      mem_type_id = H5T_NATIVE_DOUBLE; break;
    case 4:
      mem_type_id = H5T_NATIVE_LONG; break;
    case 5:
      mem_type_id = H5T_NATIVE_SHORT; break;
    default:
      throw UnrecognizedParameterException("unrecognized data type \"" + data_type + "\"");
    }
    if (dataset_id != H5I_INVALID_HID) {
      status = H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, dxpl_id, buffer);
      if (status == 0) status = H5Dclose(dataset_id);
      else status = (herr_t)-2;
    } else {
      status = (herr_t)-1;
    }
  } else {
    throw UnrecognizedParameterException("unrecognized data type \"" + data_type + "\"");
  }
  return status;
}
/* >>> End of HDFFile class implementation <<< */

/* >>> VirtualAsciiFile class implementation <<< */
VirtualAsciiFile::VirtualAsciiFile(int32_t lines) {
  _file_lines = new vector<string>();
  for (int32_t li = 0; li < lines; li++) {
    _file_lines->push_back("");
  }
}

VirtualAsciiFile::VirtualAsciiFile(const VirtualAsciiFile& rhs) {
  // _num_lines = rhs._num_lines;
  _file_lines = new vector<string>();
  for (vector<string>::iterator it = rhs._file_lines->begin(); it != rhs._file_lines->end(); ++it) {
    _file_lines->push_back(*it);
  }
}

#ifdef MPI_VERSION
VirtualAsciiFile::VirtualAsciiFile(const mixMPI *mpidata, int rr) {
  // receive _num_lines from MPI process rr
  int32_t num_lines;
  MPI_Recv(&num_lines, 1, MPI_INT32_T, rr, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  int32_t mysize;
  _file_lines = new vector<string>();
  // loop over data to receive
  for (int32_t zi=0; zi<num_lines; zi++) {
    // receive the size of the string to receive
    MPI_Recv(&mysize, 1, MPI_INT32_T, rr, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // allocate the buffer accordingly
    char buffer[mysize+1];
    // receive the char buffer
    MPI_Recv(buffer, mysize+1, MPI_CHAR, rr, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // append to the vector
    _file_lines->push_back(buffer);
  }
}
#endif

VirtualAsciiFile::~VirtualAsciiFile() {
  // is it necessary to pop them out one by one? isn't there the dedicated method of std::vector to clean the vector?
  // besides, shouldn't this be done anyway by the destructor of std:vector?
  _file_lines->clear();
  // while (_file_lines->size() > 0) {
  //   _file_lines->pop_back();
  // }
  if (_file_lines != NULL) delete _file_lines;
}

void VirtualAsciiFile::append(VirtualAsciiFile& rhs) {
  // concatenate the virtualasciifile pointed by rhs to the current one
  // can't we use the dedicated method insert of std::vector to do the appending, instead of an explicit loop?
  for (vector<string>::iterator it = rhs._file_lines->begin(); it != rhs._file_lines->end(); ++it) {
    _file_lines->push_back(*it);
  }
}

void VirtualAsciiFile::append_line(const string& line) {
  // would it be worth reimplementing a sprintf-like method, so that we can give it all the arguments we would give to sprintf and get rid of the intermediate buffer completely?
  // append a line of output to the virtualasciifile
  _file_lines->push_back(line);
}

int VirtualAsciiFile::append_to_disk(const std::string& file_name) {
  // dump to disk the contents of the virtualasciifile, appending at the end of the given file_name
  int result = 0;
  if (_file_lines->size() > 0) {
    fstream output_file;
    output_file.open(file_name, ios::app);
    if (output_file.is_open()) {
      for (vector<string>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
	output_file << *it;
      }
      output_file.close();
    } else {
      result = 1;
    }
  }
  return result;
}

int VirtualAsciiFile::insert(int32_t position, VirtualAsciiFile& rhs, int32_t start, int32_t end) {
  int result = 0;
  if (start == 0 && end == 0) {
    end = rhs.number_of_lines();
  }
  int32_t final_index = position + end - start;
  if (final_index <= number_of_lines()) {
    for (int32_t li = start; li < end; li++) {
      // since here we are replacing the previous placeholder empty strings, make sure they are properly released when they are replaced (i.e. try it with a simple hello world example and pass it through valgrind)
      _file_lines->at(position++) = rhs._file_lines->at(li);
    }
  } else {
    // ERROR: target file is too long;
    result = 1;
  }
  return result;
}

int VirtualAsciiFile::write_to_disk(const std::string& file_name) {
  // dump to disk the contents of the virtualasciifile, replacing the given file_name
  int result = 0;
  if (_file_lines->size() > 0) {
    fstream output_file;
    output_file.open(file_name, ios::out);
    if (output_file.is_open()) {
      for (vector<string>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
	output_file << *it;
      }
      output_file.close();
    } else {
      result = 1;
    }
  }
  return result;
}

#ifdef MPI_VERSION
void VirtualAsciiFile::mpisend(const mixMPI *mpidata) {
  // Send VirtualAsciiFile instance to MPI process 0 via MPISend() calls
  // first send the size
  int32_t num_lines =  _file_lines->size();
  MPI_Send(&num_lines, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
  // now loop over data to send
  if (num_lines>0) {
    int32_t mysize;
    for (vector<string>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
      // send the length of each string
      mysize = (int32_t) it->size();
      MPI_Send(&mysize, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      // send the string itself
      MPI_Send(it->c_str(), mysize+1, MPI_CHAR, 0, 10, MPI_COMM_WORLD);
    }
  }
}
#endif

/* >>> End of VirtualAsciiFile class implementation <<< */

// /* >>> VirtualBinaryLine class implementation <<< */
VirtualBinaryLine::VirtualBinaryLine(int mydata) {
  _data_size = sizeof(mydata);
  int *buffer = (int *) malloc(_data_size);
  *buffer = mydata;
  _data_pointer = reinterpret_cast<char *>(buffer);
}

VirtualBinaryLine::VirtualBinaryLine(double mydata) {
  _data_size = sizeof(mydata);
  double *buffer = (double *) malloc(_data_size);
  *buffer = mydata;
  _data_pointer = reinterpret_cast<char *>(buffer);
}

VirtualBinaryLine::VirtualBinaryLine(float mydata) {
  _data_size = sizeof(mydata);
  float *buffer = (float *) malloc(_data_size);
  *buffer = mydata;
  _data_pointer = reinterpret_cast<char *>(buffer);
}

VirtualBinaryLine::VirtualBinaryLine(long mydata) {
  _data_size = sizeof(mydata);
  long *buffer = (long *) malloc(_data_size);
  *buffer = mydata;
  _data_pointer = reinterpret_cast<char *>(buffer);
}

VirtualBinaryLine::VirtualBinaryLine(dcomplex mydata) {
  _data_size = sizeof(mydata);
  dcomplex *buffer = (dcomplex *) malloc(_data_size);
  *buffer = mydata;
  _data_pointer = reinterpret_cast<char *>(buffer);
}

// VirtualBinaryLine::VirtualBinaryLine(complex mydata) {
//   _data_size = sizeof(mydata);
//   complex *buffer = (complex *) malloc(_data_size);
//   *buffer = mydata;
//   _data_pointer = reinterpret_cast<char *>(buffer);
// }

VirtualBinaryLine::VirtualBinaryLine(const VirtualBinaryLine& rhs) {
  _data_size = rhs._data_size;
  _data_pointer = reinterpret_cast<char *>(malloc(rhs._data_size));
  memcpy(_data_pointer, rhs._data_pointer, _data_size);
}

#ifdef MPI_VERSION
VirtualBinaryLine::VirtualBinaryLine(const mixMPI *mpidata, int rr) {
  // receive mysize from MPI process rr
  int32_t mysize;
  MPI_Recv(&mysize, 1, MPI_INT32_T, rr, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  _data_size = mysize;
  // allocate the buffer accordingly
  _data_pointer = reinterpret_cast<char *>(malloc(mysize));
  // receive the char buffer
  MPI_Recv(_data_pointer, mysize, MPI_CHAR, rr, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

#endif

VirtualBinaryLine::~VirtualBinaryLine() {
  if (_data_pointer != NULL) {
    free(_data_pointer);
    _data_pointer = NULL;
  }
}

#ifdef MPI_VERSION
void VirtualBinaryLine::mpisend(const mixMPI *mpidata) {
  // Send VirtualBinaryLine instance to MPI process 0 via MPISend() calls
  // first send the size
  int32_t mysize =  _data_size;
  MPI_Send(&mysize, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
  // now send the data
  MPI_Send(_data_pointer, mysize, MPI_CHAR, 0, 10, MPI_COMM_WORLD);
}
#endif
/* >>> End of VirtualBinaryLine class implementation <<< */


/* >>> VirtualBinaryFile class implementation <<< */
VirtualBinaryFile::VirtualBinaryFile() {
  _file_lines = new vector<VirtualBinaryLine>();
}

VirtualBinaryFile::VirtualBinaryFile(const VirtualBinaryFile& rhs) {
  // _num_lines = rhs._num_lines;
  _file_lines = new vector<VirtualBinaryLine>();
  for (vector<VirtualBinaryLine>::iterator it = rhs._file_lines->begin(); it != rhs._file_lines->end(); ++it) {
    _file_lines->push_back(VirtualBinaryLine(*it));
  }
}

#ifdef MPI_VERSION
VirtualBinaryFile::VirtualBinaryFile(const mixMPI *mpidata, int rr) {
  // receive _num_lines from MPI process rr
  int32_t num_lines;
  MPI_Recv(&num_lines, 1, MPI_INT32_T, rr, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  _file_lines = new vector<VirtualBinaryLine>();
  // loop over data to receive
  for (int32_t zi=0; zi<num_lines; zi++) {
    // receive the line of data
    _file_lines->push_back(VirtualBinaryLine(mpidata, rr));
  }
}
#endif

VirtualBinaryFile::~VirtualBinaryFile() {
  // is it necessary to pop them out one by one? isn't there the dedicated method of std::vector to clean the vector?
  // besides, shouldn't this be done anyway by the destructor of std:vector?
  // for (vector<VirtualBinaryLine>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
  //   delete it;
  // }
  _file_lines->clear();
  // while (_file_lines->size() > 0) {
  //   _file_lines->pop_back();
  // }
  if (_file_lines != NULL) delete _file_lines;
}

void VirtualBinaryFile::append(VirtualBinaryFile& rhs) {
  // concatenate the virtualasciifile pointed by rhs to the current one
  // can't we use the dedicated method insert of std::vector to do the appending, instead of an explicit loop?
  for (vector<VirtualBinaryLine>::iterator it = rhs._file_lines->begin(); it != rhs._file_lines->end(); ++it) {
    _file_lines->push_back(VirtualBinaryLine(*it));
  }
}

void VirtualBinaryFile::append_line(const VirtualBinaryLine& line) {
  // would it be worth reimplementing a sprintf-like method, so that we can give it all the arguments we would give to sprintf and get rid of the intermediate buffer completely?
  // append a line of output to the virtualasciifile
  _file_lines->push_back(VirtualBinaryLine(line));
}

int VirtualBinaryFile::append_to_disk(const std::string& file_name) {
  // dump to disk the contents of the virtualasciifile, appending at the end of the given file_name
  int result = 0;
  if (_file_lines->size() > 0) {
    fstream output_file;
    output_file.open(file_name, ios::app | ios::binary);
    if (output_file.is_open()) {
      for (vector<VirtualBinaryLine>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
	output_file.write(it->_data_pointer, it->_data_size);
      }
      output_file.close();
    } else {
      result = 1;
    }
  }
  return result;
}

// int VirtualBinaryFile::insert(int32_t position, VirtualBinaryFile& rhs, int32_t start, int32_t end) {
//   int result = 0;
//   if (start == 0 && end == 0) {
//     end = rhs.number_of_lines();
//   }
//   int32_t final_index = position + end - start;
//   if (final_index <= number_of_lines()) {
//     for (int32_t li = start; li < end; li++) {
//       // since here we are replacing the previous placeholder empty strings, make sure they are properly released when they are replaced (i.e. try it with a simple hello world example and pass it through valgrind)
//       VirtualBinaryLine templine = VirtualBinaryLine(rhs._file_lines->at(li));
//       _file_lines->at(position++) = templine;
//     }
//   } else {
//     // ERROR: target file is too long;
//     result = 1;
//   }
//   return result;
// }

int VirtualBinaryFile::write_to_disk(const std::string& file_name) {
  // dump to disk the contents of the virtualasciifile, replacing the given file_name
  int result = 0;
  if (_file_lines->size() > 0) {
    fstream output_file;
    output_file.open(file_name, ios::out | ios::binary);
    if (output_file.is_open()) {
      for (vector<VirtualBinaryLine>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
	output_file.write(it->_data_pointer, it->_data_size);
      }
      output_file.close();
    } else {
      result = 1;
    }
  }
  return result;
}

#ifdef MPI_VERSION
void VirtualBinaryFile::mpisend(const mixMPI *mpidata) {
  // Send VirtualBinaryFile instance to MPI process 0 via MPISend() calls
  // first send the size
  int32_t num_lines =  _file_lines->size();
  MPI_Send(&num_lines, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
  // now loop over data to send
  if (num_lines>0) {
    for (vector<VirtualBinaryLine>::iterator it = _file_lines->begin(); it != _file_lines->end(); ++it) {
      it->mpisend(mpidata);
    }
  }
}
#endif

/* >>> End of VirtualBinaryFile class implementation <<< */
