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

/*! \file Configuration.cpp
 *
 * \brief Implementation of the configuration classes.
 */

#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <hdf5.h>
#include <regex>
#include <string>
#include <string.h>
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

#ifndef INCLUDE_PARSERS_H_
#include "../include/Parsers.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H
#include "../include/Commons.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

using namespace std;

GeometryConfiguration::GeometryConfiguration(
  int nsph, int lm, int in_pol, int npnt, int npntts, int isam,
  int li, int le, np_int mxndm, int iavm, double *x, double *y,
  double *z, double in_th_start, double in_th_step, double in_th_end,
  double sc_th_start, double sc_th_step, double sc_th_end,
  double in_ph_start, double in_ph_step, double in_ph_end,
  double sc_ph_start, double sc_ph_step, double sc_ph_end, int jwtm
) {
  _number_of_spheres = nsph;
  _l_max = lm;
  _in_pol = in_pol;
  _npnt = npnt;
  _npntts = npntts;
  _isam = isam;
  _li = li;
  _le = le;
  _mxndm = mxndm;
  _iavm = iavm;
  _in_theta_start = in_th_start;
  _in_theta_step = in_th_step;
  _in_theta_end = in_th_end;
  _in_phi_start = in_ph_start;
  _in_phi_step = in_ph_step;
  _in_phi_end = in_ph_end;
  _sc_theta_start = sc_th_start;
  _sc_theta_step = sc_th_step;
  _sc_theta_end = sc_th_end;
  _sc_phi_start = sc_ph_start;
  _sc_phi_step = sc_ph_step;
  _sc_phi_end = sc_ph_end;
  _jwtm = jwtm;
  _sph_x = x;
  _sph_y = y;
  _sph_z = z;
  _refine_flag = 0;
  _dyn_order_flag = 1;
  _host_ram_gb = 0.0;
  _gpu_ram_gb = 0.0;
}

GeometryConfiguration::GeometryConfiguration(const GeometryConfiguration& rhs)
{
  _number_of_spheres = rhs._number_of_spheres;
  _l_max = rhs._l_max;
  _in_pol = rhs._in_pol;
  _npnt = rhs._npnt;
  _npntts = rhs._npntts;
  _isam = rhs._isam;
  _li = rhs._li;
  _le = rhs._le;
  _mxndm = rhs._mxndm;
  _iavm = rhs._iavm;
  _in_theta_start = rhs._in_theta_start;
  _in_theta_step = rhs._in_theta_step;
  _in_theta_end = rhs._in_theta_end;
  _in_phi_start = rhs._in_phi_start;
  _in_phi_step = rhs._in_phi_step;
  _in_phi_end = rhs._in_phi_end;
  _sc_theta_start = rhs._sc_theta_start;
  _sc_theta_step = rhs._sc_theta_step;
  _sc_theta_end = rhs._sc_theta_end;
  _sc_phi_start = rhs._sc_phi_start;
  _sc_phi_step = rhs._sc_phi_step;
  _sc_phi_end = rhs._sc_phi_end;
  _jwtm = rhs._jwtm;
  _sph_x = new double[_number_of_spheres]();
  _sph_y = new double[_number_of_spheres]();
  _sph_z = new double[_number_of_spheres]();
  for (int ni=0; ni < _number_of_spheres; ni++) {
    _sph_x[ni] = rhs._sph_x[ni];
    _sph_y[ni] = rhs._sph_y[ni];
    _sph_z[ni] = rhs._sph_z[ni];
  }
  _refine_flag = rhs._refine_flag;
  _dyn_order_flag = rhs._dyn_order_flag;
  _host_ram_gb = rhs._host_ram_gb;
  _gpu_ram_gb = rhs._gpu_ram_gb;
}

#ifdef MPI_VERSION
GeometryConfiguration::GeometryConfiguration(const mixMPI *mpidata) {
  MPI_Bcast(&_number_of_spheres, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_l_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_pol, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npntts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_isam, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_li, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_le, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // I have to jump through some hoops because the size of np_int is not fixed a priori
  char *byte_mxndm = (char *) &_mxndm;
  MPI_Bcast(byte_mxndm, sizeof(np_int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_iavm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_theta_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_theta_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_theta_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_phi_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_phi_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_phi_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_theta_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_theta_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_theta_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_phi_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_phi_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_phi_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_jwtm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  _sph_x = new double[_number_of_spheres]();
  _sph_y = new double[_number_of_spheres]();
  _sph_z = new double[_number_of_spheres]();
  MPI_Bcast(_sph_x, _number_of_spheres, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_sph_y, _number_of_spheres, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_sph_z, _number_of_spheres, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_refine_flag, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_dyn_order_flag, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_host_ram_gb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_gpu_ram_gb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void GeometryConfiguration::mpibcast(const mixMPI *mpidata) {
  MPI_Bcast(&_number_of_spheres, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_l_max, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_pol, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npntts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_isam, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_li, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_le, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // I have to jump through some hoops because the size of np_int is not fixed a priori
  char *byte_mxndm = (char *) &_mxndm;
  MPI_Bcast(byte_mxndm, sizeof(np_int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_iavm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_theta_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_theta_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_theta_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_phi_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_phi_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_in_phi_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_theta_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_theta_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_theta_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_phi_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_phi_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_sc_phi_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_jwtm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(_sph_x, _number_of_spheres, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_sph_y, _number_of_spheres, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_sph_z, _number_of_spheres, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_refine_flag, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_dyn_order_flag, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_host_ram_gb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_gpu_ram_gb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
#endif

GeometryConfiguration::~GeometryConfiguration() {
  delete[] _sph_x;
  delete[] _sph_y;
  delete[] _sph_z;
}

GeometryConfiguration* GeometryConfiguration::from_legacy(const std::string& file_name) {
  int num_lines = 0;
  int last_read_line = 0;
  string *file_lines;
  string str_target, str_num;
  smatch m;
  try {
    file_lines = load_file(file_name, &num_lines);
  } catch (...) {
    OpenConfigurationFileException ex(file_name);
    throw ex;
  }
  // Read the legacy FORTRAN mandatory configuration data
  int _nsph = 0, _lm = 0, _in_pol = 0, _npnt = 0, _npntts = 0, _isam = 0;
  int _li = 0, _le = 0, _iavm = 0, num_params = 0;
  np_int _mxndm = 0;
  regex re = regex("-?[0-9]+");
  str_target = file_lines[last_read_line];
  while(regex_search(str_target, m, re)) {
    str_target = m.suffix().str();
    num_params++;
  }
  str_target = file_lines[last_read_line++];
  regex_search(str_target, m, re);
  _nsph = stoi(m.str());
  if (num_params == 6) {
    for (int ri = 0; ri < 5; ri++) {
      str_target = m.suffix().str();
      regex_search(str_target, m, re);
      if (ri == 0) _lm = stoi(m.str());
      if (ri == 1) _in_pol = stoi(m.str());
      if (ri == 2) _npnt = stoi(m.str());
      if (ri == 3) _npntts = stoi(m.str());
      if (ri == 4) _isam = stoi(m.str());
    }
  } else if (num_params == 9) {
    for (int ri = 0; ri < 8; ri++) {
      str_target = m.suffix().str();
      regex_search(str_target, m, re);
      if (ri == 0) _li = stoi(m.str());
      if (ri == 1) _le = stoi(m.str());
      if (ri == 2) _mxndm = stoi(m.str());
      if (ri == 3) _in_pol = stoi(m.str());
      if (ri == 4) _npnt = stoi(m.str());
      if (ri == 5) _npntts = stoi(m.str());
      if (ri == 6) _iavm = stoi(m.str());
      if (ri == 7) _isam = stoi(m.str());
    }
  } else {
    OpenConfigurationFileException ex("ERROR: " + file_name + " is not a recognized input file.");
    throw ex;
  }
  double *x, *y, *z;
  x = new double[_nsph];
  y = new double[_nsph];
  z = new double[_nsph];
  if (_nsph == 1) {
    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
  } else {
    for (int i = 0; i < _nsph; i++) {
      str_target = file_lines[last_read_line++];
      re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
      for (int ri = 0; ri < 3; ri++) {
	regex_search(str_target, m, re);
	str_num = regex_replace(m.str(), regex("D"), "e");
	str_num = regex_replace(str_num, regex("d"), "e");
	if (ri == 0) x[i] = stod(str_num);
	if (ri == 1) y[i] = stod(str_num);
	if (ri == 2) z[i] = stod(str_num);
	str_target = m.suffix().str();
      }
    }
  }
  double in_th_start, in_th_end, in_th_step, sc_th_start, sc_th_end, sc_th_step;
  re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
  str_target = file_lines[last_read_line++];
  for (int ri = 0; ri < 6; ri++) {
    regex_search(str_target, m, re);
    str_num = regex_replace(m.str(), regex("D"), "e");
    str_num = regex_replace(str_num, regex("d"), "e");
    if (ri == 0) in_th_start = stod(str_num);
    if (ri == 1) in_th_step = stod(str_num);
    if (ri == 2) in_th_end = stod(str_num);
    if (ri == 3) sc_th_start = stod(str_num);
    if (ri == 4) sc_th_step = stod(str_num);
    if (ri == 5) sc_th_end = stod(str_num);
    str_target = m.suffix().str();
  }

  double in_ph_start, in_ph_end, in_ph_step, sc_ph_start, sc_ph_end, sc_ph_step;
  str_target = file_lines[last_read_line++];
  for (int ri = 0; ri < 6; ri++) {
    regex_search(str_target, m, re);
    str_num = regex_replace(m.str(), regex("D"), "e");
    str_num = regex_replace(str_num, regex("d"), "e");
    if (ri == 0) in_ph_start = stod(str_num);
    if (ri == 1) in_ph_step = stod(str_num);
    if (ri == 2) in_ph_end = stod(str_num);
    if (ri == 3) sc_ph_start = stod(str_num);
    if (ri == 4) sc_ph_step = stod(str_num);
    if (ri == 5) sc_ph_end = stod(str_num);
    str_target = m.suffix().str();
  }

  int fjwtm;
  re = regex("[0-9]+");
  str_target = file_lines[last_read_line++];
  regex_search(str_target, m, re);
  fjwtm = stoi(m.str());
  // Mandatory configuration data were read. Create the configuration object.
  GeometryConfiguration *conf = new GeometryConfiguration(
    _nsph, _lm, _in_pol, _npnt, _npntts, _isam,
    _li, _le, _mxndm, _iavm, x, y, z,
    in_th_start, in_th_step, in_th_end,
    sc_th_start, sc_th_step, sc_th_end,
    in_ph_start, in_ph_step, in_ph_end,
    sc_ph_start, sc_ph_step, sc_ph_end,
    fjwtm
  );

  // Read optional configuration data used only by the C++ code.
  while (num_lines > last_read_line) {
    str_target = file_lines[last_read_line++];
    if (str_target.size() > 15) {
      if (str_target.substr(0, 15).compare("USE_REFINEMENT=") == 0) {
	re = regex("[0-9]+");
	regex_search(str_target, m, re);
	short refine_flag = (short)stoi(m.str());
	conf->_refine_flag = refine_flag;
      }
      if (str_target.substr(0, 15).compare("USE_DYN_ORDERS=") == 0) {
	re = regex("[0-9]+");
	regex_search(str_target, m, re);
	short dyn_order_flag = (short)stoi(m.str());
	conf->_dyn_order_flag = dyn_order_flag;
      }
    }
    if (str_target.size() > 12) {
      if (str_target.substr(0, 12).compare("HOST_RAM_GB=") == 0) {
	double ram_gb = (double)stod(str_target.substr(12, str_target.length()));
	conf->_host_ram_gb = ram_gb;
      }
    }
    if (str_target.size() > 11) {
      if (str_target.substr(0, 11).compare("GPU_RAM_GB=") == 0) {
	double ram_gb = (double)stod(str_target.substr(11, str_target.length()));
	conf->_gpu_ram_gb = ram_gb;
      }
    }
  }
  
  // Clean up memory and return configuration object.
  delete[] file_lines;
  return conf;
}

ScattererConfiguration::ScattererConfiguration(
  int nsph, int configs, double *scale_vector, int nxi,
  const std::string& variable_name, int *iog_vector,
  double *ros_vector, int *nshl_vector, double **rcf_vector,
  int dielectric_func_type, dcomplex ***dc_matrix, bool is_external,
  double ex, double w, double x
) {
  _number_of_spheres = nsph;
  _configurations = configs;
  _number_of_scales = nxi;
  _reference_variable_name = variable_name;
  _iog_vec = iog_vector;
  _radii_of_spheres = ros_vector;
  _nshl_vec = nshl_vector;
  _use_external_sphere = is_external;
  _max_layers = 0;
  for (int li = 0; li < _configurations; li++) {
    if (_max_layers < _nshl_vec[li]) {
      _max_layers = _nshl_vec[li];
      if (li == 0 && _use_external_sphere) _max_layers += 1;
    }
  }
  _rcf = rcf_vector;
  _idfc = dielectric_func_type;
  _dc0_matrix = dc_matrix;
  _exdc = ex;
  _wp = w;
  _xip = x;
  _scale_vec = new double[_number_of_scales]();
  if (variable_name == "XIV") {
    for (int xi = 0; xi < nxi; xi++) _scale_vec[xi] = scale_vector[xi];
  } else {
    const double pi2 = 2.0 * acos(-1.0);
    const double evc = 6.5821188e-16;
    for (int si = 0; si < number_of_scales; si++) {
      double value = scale_vector[si];
      if (variable_name.compare("WNS") == 0) value *= (3.0e8 / wp);
      else if (variable_name.compare("WLS") == 0) value = pi2 / value * 3.0e8 / wp;
      else if (variable_name.compare("PUS") == 0) value /= wp;
      else if (variable_name.compare("EVS") == 0) value /= (evc * wp);
      _scale_vec[si] = value;
    }
  }
}
  
ScattererConfiguration::ScattererConfiguration(const ScattererConfiguration& rhs)
{
  _number_of_spheres = rhs._number_of_spheres;
  _configurations = rhs._configurations;
  _number_of_scales = rhs._number_of_scales;
  _reference_variable_name = rhs._reference_variable_name;
  _idfc = rhs._idfc; 
  _use_external_sphere = rhs._use_external_sphere;
  _exdc = rhs._exdc;
  _wp = rhs._wp;
  _xip = rhs._xip;
  _max_layers = rhs._max_layers;
  _iog_vec = new int[_number_of_spheres]();
  _radii_of_spheres = new double[_number_of_spheres]();
  _nshl_vec = new int[_configurations]();
  _rcf = new double*[_configurations];
  _scale_vec = new double[_number_of_scales]();
  _dc0_matrix = new dcomplex**[_max_layers];
  for (int si = 0; si < _number_of_scales; si++) _scale_vec[si] = rhs._scale_vec[si];
  for (int si = 0; si < _number_of_spheres; si++) {
    _iog_vec[si] = rhs._iog_vec[si];
  }
  int dim3 = (_idfc == 0) ? _number_of_scales : 1;
  for (int si = 0; si < _configurations; si++) {
    _radii_of_spheres[si] = rhs._radii_of_spheres[si];
    _nshl_vec[si] = rhs._nshl_vec[si];
    _rcf[si] = new double[_nshl_vec[si]]();
    for (int sj = 0; sj < _nshl_vec[si]; sj++) _rcf[si][sj] = rhs._rcf[si][sj];
  }
  for (int li = 0; li < _max_layers; li++) {
    _dc0_matrix[li] = new dcomplex*[_number_of_spheres];
    for (int lj = 0; lj < _number_of_spheres; lj++) {
      _dc0_matrix[li][lj] = new dcomplex[dim3]();
      for (int lk = 0; lk < dim3; lk++) _dc0_matrix[li][lj][lk] = rhs._dc0_matrix[li][lj][lk];
    }
  }
}

#ifdef MPI_VERSION
ScattererConfiguration::ScattererConfiguration(const mixMPI *mpidata)
{
  MPI_Bcast(&_number_of_spheres, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_configurations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_number_of_scales, 1, MPI_INT, 0, MPI_COMM_WORLD);
  int itemp;
  MPI_Bcast(&itemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  char *ctemp = new char[itemp];
  MPI_Bcast(ctemp, itemp, MPI_CHAR, 0, MPI_COMM_WORLD);
  _reference_variable_name = ctemp;
  delete[] ctemp;
  MPI_Bcast(&_idfc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  itemp = sizeof(bool);
  char *ptemp = (char *) &_use_external_sphere;
  MPI_Bcast(ptemp, itemp, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_max_layers, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_exdc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_wp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_xip, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  _iog_vec = new int[_number_of_spheres]();
  MPI_Bcast(_iog_vec, _number_of_spheres, MPI_INT, 0, MPI_COMM_WORLD);
  _radii_of_spheres = new double[_configurations]();
  MPI_Bcast(_radii_of_spheres, _configurations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  _nshl_vec = new int[_configurations]();
  MPI_Bcast(_nshl_vec, _configurations, MPI_INT, 0, MPI_COMM_WORLD);
  _rcf = new double*[_configurations];
  _scale_vec = new double[_number_of_scales]();
  _dc0_matrix = new dcomplex**[_max_layers];
  MPI_Bcast(_scale_vec, _number_of_scales, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  int dim3 = (_idfc == 0) ? _number_of_scales : 1;
  for (int si = 0; si < _configurations; si++) {
    _rcf[si] = new double[_nshl_vec[si]]();
    MPI_Bcast(_rcf[si], _nshl_vec[si], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  for (int li = 0; li < _max_layers; li++) {
    _dc0_matrix[li] = new dcomplex*[_number_of_spheres];
    for (int lj = 0; lj < _number_of_spheres; lj++) {
      _dc0_matrix[li][lj] = new dcomplex[dim3]();
      MPI_Bcast(_dc0_matrix[li][lj], dim3, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    }
  }
}

void ScattererConfiguration::mpibcast(const mixMPI *mpidata) {
  MPI_Bcast(&_number_of_spheres, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_configurations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_number_of_scales, 1, MPI_INT, 0, MPI_COMM_WORLD);
  int itemp = _reference_variable_name.length()+1;
  MPI_Bcast(&itemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  char *ctemp = strdup(_reference_variable_name.c_str());
  MPI_Bcast(ctemp, itemp, MPI_CHAR, 0, MPI_COMM_WORLD);
  free(ctemp);
  MPI_Bcast(&_idfc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  itemp = sizeof(bool);
  char *ptemp = (char *) &_use_external_sphere;
  MPI_Bcast(ptemp, itemp, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_max_layers, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_exdc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_wp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_xip, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_iog_vec, _number_of_spheres, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(_radii_of_spheres, _configurations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(_nshl_vec, _configurations, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(_scale_vec, _number_of_scales, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  int dim3 = (_idfc == 0) ? _number_of_scales : 1;
  for (int si = 0; si < _configurations; si++) {
    MPI_Bcast(_rcf[si], _nshl_vec[si], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  for (int li = 0; li < _max_layers; li++) {
    for (int lj = 0; lj < _number_of_spheres; lj++) {
      MPI_Bcast(_dc0_matrix[li][lj], dim3, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    }
  }
}
#endif

ScattererConfiguration::~ScattererConfiguration() {
  for (int i = 0; i < _max_layers; i++) {
    for (int j = 0; j < _number_of_spheres; j++) {
      delete[] _dc0_matrix[i][j];
    }
    delete[] _dc0_matrix[i];
  }
  delete[] _dc0_matrix;
  for (int i = 0; i < _configurations; i++) {
    delete[] _rcf[i];
  }
  delete[] _rcf;
  delete[] _nshl_vec;
  delete[] _radii_of_spheres;
  delete[] _iog_vec;
  delete[] _scale_vec;
}

ScattererConfiguration* ScattererConfiguration::from_binary(const std::string& file_name, const std::string& mode) {
  ScattererConfiguration *conf = NULL;
  if (mode.compare("LEGACY") == 0) {
    conf = ScattererConfiguration::from_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    conf = ScattererConfiguration::from_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"";
    throw UnrecognizedConfigurationException(message);
  }
  return conf;
}

ScattererConfiguration* ScattererConfiguration::from_dedfb(const std::string& dedfb_file_name) {
  int num_lines = 0;
  int last_read_line = 0;
  int last_configuration;
  string *file_lines;
  regex re;
  smatch m;
  try {
    file_lines = load_file(dedfb_file_name, &num_lines);
  } catch (...) {
    OpenConfigurationFileException ex(dedfb_file_name);
    throw ex;
  }
  int fnsph, fies;
  re = regex("[0-9]+");
  string str_target = file_lines[last_read_line];
  for (int ri = 0; ri < 2; ri++) {
    regex_search(str_target, m, re);
    if (ri == 0) fnsph = stoi(m.str());
    if (ri == 1) fies = stoi(m.str());
    str_target = m.suffix().str();
  }
  if (fies != 0) fies = 1;
  double fexdc, fwp, fxip;
  str_target = file_lines[++last_read_line];
  re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
  for (int ri = 0; ri < 3; ri++) {
    regex_search(str_target, m, re);
    string str_number = m.str();
    str_number = regex_replace(str_number, regex("D"), "e");
    str_number = regex_replace(str_number, regex("d"), "e");
    if (ri == 0) fexdc = stod(str_number);
    if (ri == 1) fwp = stod(str_number);
    if (ri == 2) fxip = stod(str_number);
    str_target = m.suffix().str();
  }
  int fidfc, fnxi, finstpc, finsn;
  re = regex("-?[0-9]+");
  for (int ri = 0; ri < 4; ri++) {
    regex_search(str_target, m, re);
    if (ri == 0) fidfc = stoi(m.str());
    if (ri == 1) fnxi = stoi(m.str());
    if (ri == 2) finstpc = stoi(m.str());
    if (ri == 3) finsn = stoi(m.str());
    str_target = m.suffix().str();
  }

  double *variable_vector;
  string variable_name;

  if (fidfc < 0) { // Diel. functions at XIP value and XI is scale factor
    variable_name = "XIV";
    if (finstpc < 1) { // The variable vector is explicitly defined.
      double xi;
      List<double> *xi_vector = new List<double>(1);
      str_target = file_lines[++last_read_line];
      re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
      regex_search(str_target, m, re);
      string str_number = m.str();
      str_number = regex_replace(str_number, regex("D"), "e");
      str_number = regex_replace(str_number, regex("d"), "e");
      xi = stod(str_number);
      xi_vector->set(0, xi);
      for (int jxi310 = 1; jxi310 < fnxi; jxi310++) {
	str_target = file_lines[++last_read_line];
	regex_search(str_target, m, re);
	str_number = m.str();
	str_number = regex_replace(str_number, regex("D"), "e");
	str_number = regex_replace(str_number, regex("d"), "e");
	xi = stod(str_number);
	xi_vector->append(xi);
      }
      variable_vector = xi_vector->to_array();
      delete xi_vector;
    } else { // instpc >= 1: the variable vector is defined in steps
      double xi, xi_step;
      str_target = file_lines[++last_read_line];
      re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
      regex_search(str_target, m, re);
      for (int ri = 0; ri < 2; ri++) {
	regex_search(str_target, m, re);
	string str_number = m.str();
	str_number = regex_replace(str_number, regex("D"), "e");
	str_number = regex_replace(str_number, regex("d"), "e");
	if (ri == 0) xi = stod(str_number);
	if (ri == 1) xi_step = stod(str_number);
	str_target = m.suffix().str();
      }
      variable_vector = new double[fnxi]();
      for (int jxi320 = 0; jxi320 < fnxi; jxi320++) {
	variable_vector[jxi320] = xi;
	xi += xi_step;
      }
    }
  } else { // idfc >= 0
    variable_vector = new double[fnxi]();
    if (finstpc == 0) { // The variable vector is explicitly defined
      double vs;
      for (int jxi_r = 0; jxi_r < fnxi; jxi_r++) {
	str_target = file_lines[++last_read_line];
	re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
	regex_search(str_target, m, re);
	string str_number = m.str();
	str_number = regex_replace(str_number, regex("D"), "e");
	str_number = regex_replace(str_number, regex("d"), "e");
	vs = stod(str_number);
	variable_vector[jxi_r] = vs;
      }
      switch (finsn) {
      case 1: //xi vector definition
	variable_name = "XIV";
	break;
      case 2: //wave number vector definition
	variable_name = "WNS";
	break;
      case 3: //wavelength vector definition
	variable_name = "WLS";
	break;
      case 4: //pu vector definition
	variable_name = "PUS";
	break;
      case 5: //eV vector definition
	variable_name = "EVS";
	break;
      }
    } else { // The variable vector needs to be computed in steps
      double vs, vs_step;
      str_target = file_lines[++last_read_line];
      re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
      regex_search(str_target, m, re);
      for (int ri = 0; ri < 2; ri++) {
	regex_search(str_target, m, re);
	string str_number = m.str();
	str_number = regex_replace(str_number, regex("D"), "e");
	str_number = regex_replace(str_number, regex("d"), "e");
	if (ri == 0) vs = stod(str_number);
	if (ri == 1) vs_step = stod(str_number);
	str_target = m.suffix().str();
      }
      for (int jxi110w = 0; jxi110w < fnxi; jxi110w++) {
	variable_vector[jxi110w] = vs;
	vs += vs_step;
      }
      switch (finsn) {
      case 1: //xi vector definition
	variable_name = "XIV";
	break;
      case 2: //wave number vector definition
	variable_name = "WNS";
	break;
      case 3: //wavelength vector definition
	variable_name = "WLS";
	break;
      case 4: //pu vector definition
	variable_name = "PUS";
	break;
      case 5: //eV vector definition
	variable_name = "EVS";
	break;
      }
    }
  }
  int *iog_vector = new int[fnsph]();
  str_target = file_lines[++last_read_line];
  re = regex("[0-9]+");
  int configurations = 0;
  for (int i = 1; i <= fnsph; i++) {
    bool success = regex_search(str_target, m, re);
    if (success) {
      iog_vector[i - 1] = stoi(m.str());
      str_target = m.suffix().str();
      if (iog_vector[i - 1] >= i) configurations++;
    } else {
      str_target = file_lines[++last_read_line];
      i--;
    }
  }
  double *ros_vector = new double[configurations]();
  int *nshl_vector = new int[configurations]();
  double **rcf_vector = new double*[configurations];
  int max_layers = 0;
  last_configuration = 0;
  for (int i113 = 1; i113 <= fnsph; i113++) {
    if (iog_vector[i113 - 1] < i113) continue;
    else last_configuration++;
    str_target = file_lines[++last_read_line];
    re = regex("[0-9]+");
    regex_search(str_target, m, re);
    nshl_vector[last_configuration - 1] = stoi(m.str());
    if (max_layers < nshl_vector[last_configuration - 1])
      max_layers = nshl_vector[last_configuration - 1];
    if (last_configuration == 1 && fies == 1) max_layers += 1;
    str_target = m.suffix().str();
    re = regex("-?[0-9]+\\.[0-9]+([eEdD][-+]?[0-9]+)?");
    regex_search(str_target, m, re);
    string str_number = m.str();
    str_number = regex_replace(str_number, regex("D"), "e");
    str_number = regex_replace(str_number, regex("d"), "e");
    ros_vector[last_configuration - 1] = stod(str_number);
    int nsh = nshl_vector[last_configuration - 1];
    if (i113 == 1) nsh += fies;
    rcf_vector[last_configuration - 1] = new double[nsh]();
    for (int ns112 = 0; ns112 < nsh; ns112++) {
      str_target = file_lines[++last_read_line];
      regex_search(str_target, m, re);
      str_number = m.str();
      str_number = regex_replace(str_number, regex("D"), "e");
      str_number = regex_replace(str_number, regex("d"), "e");
      rcf_vector[last_configuration - 1][ns112] = stod(str_number);
    } // ns112 loop
  } // i113 loop

  dcomplex ***dc0m = new dcomplex**[max_layers];
  dcomplex *dc0 = new dcomplex[max_layers]();
  int dim3 = (fidfc == 0) ? fnxi : 1;
  for (int di = 0; di < max_layers; di++) {
    dc0m[di] = new dcomplex*[fnsph];
    for (int dj = 0; dj < fnsph; dj++) dc0m[di][dj] = new dcomplex[dim3]();
  } // di loop
  for (int jxi468 = 1; jxi468 <= fnxi; jxi468++) {
    if (fidfc != 0 && jxi468 > 1) continue; // jxi468 loop
    last_configuration = 0;
    for (int i162 = 1; i162 <= fnsph; i162++) {
      if (iog_vector[i162 - 1] >= i162) {
	last_configuration++;
	int nsh = nshl_vector[last_configuration - 1];
	int ici = (nsh + 1) / 2;
	if (i162 == 1) ici += fies;
	for (int ic157 = 0; ic157 < ici; ic157++) {
	  str_target = file_lines[++last_read_line];
	  regex_search(str_target, m, re);
	  string str_number = m.str();
	  str_number = regex_replace(str_number, regex("D"), "e");
	  str_number = regex_replace(str_number, regex("d"), "e");
	  double rval = stod(str_number);
	  str_target = m.suffix().str();
	  regex_search(str_target, m, re);
	  str_number = m.str();
	  str_number = regex_replace(str_number, regex("D"), "e");
	  str_number = regex_replace(str_number, regex("d"), "e");
	  double ival = stod(str_number);
	  dc0[ic157] = rval + ival * I;
	  dc0m[ic157][i162 - 1][jxi468 - 1] = dc0[ic157];
	} // ic157 loop
      }
    } // i162 loop
  } // jxi468 loop
  delete[] dc0;

  ScattererConfiguration *config = new ScattererConfiguration(
    fnsph, configurations, variable_vector, fnxi, variable_name,
    iog_vector, ros_vector, nshl_vector, rcf_vector, fidfc,
    dc0m, (fies > 0), fexdc, fwp, fxip
  );
  delete[] file_lines;
  delete[] variable_vector;
  return config;
}

ScattererConfiguration* ScattererConfiguration::from_hdf5(const std::string& file_name) {
  ScattererConfiguration *conf = NULL;
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(file_name, flags);
  herr_t status = hdf_file->get_status();
  string str_name, str_type;
  if (status == 0) {
    int nsph, ies;
    int last_configuration;
    int *iog;
    double _exdc, _wp, _xip;
    int _idfc, nxi;
    int *nshl_vector;
    double *xi_vec;
    double *ros_vector;
    double **rcf_vector;
    dcomplex ***dc0m;
    status = hdf_file->read("NSPH", "INT32_(1)", &nsph);
    status = hdf_file->read("IES", "INT32_(1)", &ies);
    status = hdf_file->read("EXDC", "FLOAT64_(1)", &_exdc);
    status = hdf_file->read("WP", "FLOAT64_(1)", &_wp);
    status = hdf_file->read("XIP", "FLOAT64_(1)", &_xip);
    status = hdf_file->read("IDFC", "INT32_(1)", &_idfc);
    status = hdf_file->read("NXI", "INT32_(1)", &nxi);
    iog = new int[nsph];
    string str_type = "INT32_(" + to_string(nsph) + ")";
    status = hdf_file->read("IOGVEC", str_type, iog);
    int configuration_count = 0;
    for (int ci = 0; ci < nsph; ci++) {
      if (iog[ci] < ci + 1) continue;
      configuration_count++;
    }
    nshl_vector = new int[configuration_count]();
    ros_vector = new double[configuration_count]();
    rcf_vector = new double*[configuration_count];
    last_configuration = 0;
    int max_layers = 0;
    for (int i115 = 1; i115 <= configuration_count; i115++) {
      if (iog[i115 - 1] < i115) continue;
      else last_configuration++;
      str_name = "NSHL_" + to_string(last_configuration);
      str_type = "INT32_(1)";
      status = hdf_file->read(str_name, str_type, (nshl_vector + last_configuration - 1));
      if (max_layers < nshl_vector[last_configuration - 1]) {
	max_layers = nshl_vector[last_configuration - 1];
	if (last_configuration == 1 && ies == 1) max_layers += 1;
      }
      str_name = "ROS_" + to_string(last_configuration);
      str_type = "FLOAT64_(1)";
      status = hdf_file->read(str_name, str_type, (ros_vector + last_configuration - 1));
      int nsh = nshl_vector[last_configuration - 1];
      rcf_vector[last_configuration - 1] = new double[nsh]();
      str_name = "RCF_" +  to_string(last_configuration);
      str_type = "FLOAT64_(" + to_string(nsh) + ")";
      status = hdf_file->read(str_name, str_type, (rcf_vector[last_configuration - 1]));
    }
    xi_vec = new double[nxi]();
    str_name = "XIVEC";
    str_type = "FLOAT64_(" + to_string(nxi) + ")";
    status = hdf_file->read(str_name, str_type, xi_vec);
      
    int dim3 = (_idfc == 0) ? nxi : 1;
    int element_size = 2 * dim3 * nsph * max_layers;
    double *elements = new double[element_size]();
    str_name = "DC0M";
    str_type = "FLOAT64_(" + to_string(element_size) + ")";
    status = hdf_file->read(str_name, str_type, elements);
    dc0m = new dcomplex**[max_layers];
    int dc_index = 0;
    for (int di = 0; di < max_layers; di++) {
      dc0m[di] = new dcomplex*[nsph];
      for (int dj = 0; dj < nsph; dj++) {
	dc0m[di][dj] = new dcomplex[dim3]();
	for (int dk = 0; dk < dim3; dk++) {
	  double rval = elements[2 * dc_index];
	  double ival = elements[2 * dc_index + 1];
	  dc_index++;
	  dc0m[di][dj][dk] = rval + ival * I;
	}
      } // dj loop
    } // di loop
    delete[] elements;
    status = hdf_file->close();
    delete hdf_file;
    conf = new ScattererConfiguration(
      nsph, configuration_count, xi_vec, nxi, "XIV",
      iog, ros_vector, nshl_vector, rcf_vector, _idfc,
      dc0m, (ies == 1), _exdc, _wp, _xip
    );
    delete[] xi_vec;
  }
  return conf;
}

ScattererConfiguration* ScattererConfiguration::from_legacy(const std::string& file_name) {
  int _nsph;
  int *_iog_vec;
  int _ies;
  int last_configuration;
  double _exdc, _wp, _xip;
  int _idfc, _nxi;
  int *_nshl_vec;
  double *_xi_vec;
  double *_ros_vec;
  double **_rcf_vec;
  dcomplex ***_dc0m;

  fstream input;
  input.open(file_name.c_str(), ios::in | ios::binary);
  input.read(reinterpret_cast<char *>(&_nsph), sizeof(int));
  input.read(reinterpret_cast<char *>(&_ies), sizeof(int));
  _iog_vec = new int[_nsph]();
  for (int i = 0; i < _nsph; i++)
    input.read(reinterpret_cast<char *>(&(_iog_vec[i])), sizeof(int));
  input.read(reinterpret_cast<char *>(&_exdc), sizeof(double));
  input.read(reinterpret_cast<char *>(&_wp), sizeof(double));
  input.read(reinterpret_cast<char *>(&_xip), sizeof(double));
  input.read(reinterpret_cast<char *>(&_idfc), sizeof(int));
  input.read(reinterpret_cast<char *>(&_nxi), sizeof(int));
  _xi_vec = new double[_nxi]();
  for (int i = 0; i < _nxi; i++)
    input.read(reinterpret_cast<char *>(&(_xi_vec[i])), sizeof(double));
  int configurations = 0;
  for (int i = 1; i <= _nsph; i++) {
    if (_iog_vec[i - 1] >= i) configurations++;
  }
  _nshl_vec = new int[configurations]();
  _ros_vec = new double[configurations]();
  _rcf_vec = new double*[configurations];
  last_configuration = 0;
  for (int i115 = 1; i115 <= _nsph; i115++) {
    if (_iog_vec[i115 - 1] < i115) continue;
    else last_configuration++;
    input.read(reinterpret_cast<char *>(&(_nshl_vec[last_configuration - 1])), sizeof(int));
    input.read(reinterpret_cast<char *>(&(_ros_vec[last_configuration - 1])), sizeof(double));
    int nsh = _nshl_vec[last_configuration - 1];
    if (i115 == 1) nsh += _ies;
    _rcf_vec[last_configuration - 1] = new double[nsh]();
    for (int nsi = 0; nsi < nsh; nsi++)
      input.read(reinterpret_cast<char *>(&(_rcf_vec[last_configuration - 1][nsi])), sizeof(double));
  }
  int max_layers = 0;
  for (int li = 0; li < configurations; li++) {
    if (max_layers < _nshl_vec[li]) max_layers = _nshl_vec[li];
    if (li == 0 && _ies == 1) max_layers += 1;
  }
  _dc0m = new dcomplex**[max_layers];
  int dim3 = (_idfc == 0) ? _nxi : 1;
  for (int di = 0; di < max_layers; di++) {
    _dc0m[di] = new dcomplex*[_nsph];
    for (int dj = 0; dj < _nsph; dj++) _dc0m[di][dj] = new dcomplex[dim3]();
  } // di loop
  for (int jxi468 = 1; jxi468 <= _nxi; jxi468++) {
    if (_idfc != 0 && jxi468 > 1) continue;
    last_configuration = 0;
    for (int i162 = 1; i162 <= _nsph; i162++) {
      if (_iog_vec[i162 - 1] < i162) continue;
      else last_configuration++;
      int nsh = _nshl_vec[last_configuration - 1];
      int ici = (nsh + 1) / 2;
      if (i162 == 1) ici = ici + _ies;
      for (int i157 = 0; i157 < ici; i157++) {
	double rval, ival;
	input.read(reinterpret_cast<char *>(&rval), sizeof(double));
	input.read(reinterpret_cast<char *>(&ival), sizeof(double));
	_dc0m[i157][i162 - 1][jxi468 - 1] = rval + ival * I;
      }
    }
  }
  input.close();
    
  ScattererConfiguration *conf = new ScattererConfiguration(
    _nsph, configurations, _xi_vec, _nxi, "XIV", _iog_vec, _ros_vec,
    _nshl_vec, _rcf_vec, _idfc, _dc0m, (_ies == 1), _exdc, _wp, _xip
  );
  delete[] _xi_vec;
  return conf;
}

double ScattererConfiguration::get_max_radius() {
  double result = 0.0;
  for (int ci = 0; ci < _configurations; ci++) {
    if (_radii_of_spheres[ci] > result)
      result = _radii_of_spheres[ci];
  }
  return result;
}

double ScattererConfiguration::get_particle_radius(GeometryConfiguration *gc) {
  double result = 0.0;
  if (_use_external_sphere) {
    result = _radii_of_spheres[0] * _rcf[0][_nshl_vec[0] - 1];
  } else {
    double dist2, max_dist;
    double max_dist2 = 0.0;
    double avgX = 0.0, avgY = 0.0, avgZ = 0.0;
#pragma omp parallel for reduction(+: avgX, avgY, avgZ)
    for (int si = 0; si < _number_of_spheres; si++) {
      avgX += gc->get_sph_x(si);
      avgY += gc->get_sph_y(si);
      avgZ += gc->get_sph_z(si);
    }
    avgX /= _number_of_spheres;
    avgY /= _number_of_spheres;
    avgZ /= _number_of_spheres;
    int configuration = 0;
    for (int si = 0; si < _number_of_spheres; si++) {
      int iogi = _iog_vec[si];
      if (iogi >= si + 1) configuration++;
      double dX = gc->get_sph_x(si) - avgX;
      double dY = gc->get_sph_y(si) - avgY;
      double dZ = gc->get_sph_z(si) - avgZ;
      double radius = _radii_of_spheres[configuration - 1];
      dist2 = dX * dX + dY * dY + dZ * dZ + radius * radius;
      if (dist2 > max_dist2) {
	max_dist2 = dist2;
      }
    }
    max_dist = sqrt(max_dist2);
    result = max_dist;
  }
  return result;
}

void ScattererConfiguration::print() {
  int ies = (_use_external_sphere)? 1 : 0;
  printf("### CONFIGURATION DATA ###\n");
  printf("NSPH  = %d\n", _number_of_spheres);
  printf("ROS   = [");
  for (int i = 0; i < _configurations; i++) printf("\t%lg", _radii_of_spheres[i]);
  printf("\t]\n");
  printf("IOG   = [");
  for (int i = 0; i < _number_of_spheres; i++) printf("\t%d", _iog_vec[i]);
  printf("\t]\n");
  printf("NSHL  = [");
  for (int i = 0; i < _configurations; i++) printf("\t%d", _nshl_vec[i]);
  printf("\t]\n");
  printf("RCF   = [\n");
  for (int i = 1; i <= _configurations; i++) {
    int nsh = _nshl_vec[i - 1];
    if (i == 1) nsh += ies;
    printf("         [");
    for (int ns = 0; ns < nsh; ns++) {
      printf("\t%lg", _rcf[i - 1][ns]);
    }
    printf("\t]\n");
  }
  printf("        ]\n");
  printf("SCALE = %s\n", _reference_variable_name.c_str());
  printf("NXI   = %d\n", _number_of_scales);
  printf("VEC   = [");
  for (int i = 0; i < _number_of_scales; i++) printf("\t%lg", _scale_vec[i]);
  printf("\t]\n");
  printf("DC0M  = [\n");
  for (int i = 0; i < _max_layers; i++) {
    printf("         [\n");
    for (int j = 0; j < _number_of_spheres; j++) {
      printf("          [");
      for (int k = 0; k < _number_of_scales; k++) {
	if (idfc != 0 and k > 0) continue;
	printf("\t%lg + i(%lg)", real(_dc0_matrix[i][j][k]), imag(_dc0_matrix[i][j][k]));
      }
      printf("\t]\n");
    }
    printf("         ]\n");
  }
  printf("        ]\n");
}

void ScattererConfiguration::write_binary(const std::string& file_name, const std::string& mode) {
  if (mode.compare("LEGACY") == 0) {
    write_legacy(file_name);
  } else if (mode.compare("HDF5") == 0) {
    write_hdf5(file_name);
  } else {
    string message = "Unknown format mode: \"" + mode + "\"";
    throw UnrecognizedConfigurationException(message);
  }
}

void ScattererConfiguration::write_hdf5(const std::string& file_name) {
  int ies = (use_external_sphere)? 1 : 0;
  int last_configuration;
  List<string> *rec_name_list = new List<string>(1);
  List<string> *rec_type_list = new List<string>(1);
  List<void *> *rec_ptr_list = new List<void *>(1);
  string str_type, str_name;
  rec_name_list->set(0, "NSPH");
  rec_type_list->set(0, "INT32_(1)");
  rec_ptr_list->set(0, &_number_of_spheres);
  rec_name_list->append("IES");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&ies);
  rec_name_list->append("IOGVEC");
  str_type = "INT32_(" + to_string(_number_of_spheres) + ")";
  rec_type_list->append(str_type);
  rec_ptr_list->append(_iog_vec);
  rec_name_list->append("EXDC");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&_exdc);
  rec_name_list->append("WP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&_wp);
  rec_name_list->append("XIP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&_xip);
  rec_name_list->append("IDFC");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&_idfc);
  rec_name_list->append("NXI");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&_number_of_scales);
  rec_name_list->append("XIVEC");
  str_type = "FLOAT64_(" + to_string(number_of_scales) + ")";
  rec_type_list->append(str_type);
  rec_ptr_list->append(_scale_vec);
  last_configuration = 0;
  for (int i115 = 1; i115 <= _number_of_spheres; i115++) {
    if (_iog_vec[i115 - 1] < i115) continue;
    else last_configuration++;
    str_name = "NSHL_" + to_string(last_configuration);
    rec_name_list->append(str_name);
    rec_type_list->append("INT32_(1)");
    rec_ptr_list->append(&(_nshl_vec[last_configuration - 1])); // was not from IOG
    str_name = "ROS_" + to_string(last_configuration);
    rec_name_list->append(str_name);
    rec_type_list->append("FLOAT64_(1)");
    rec_ptr_list->append(&(_radii_of_spheres[last_configuration - 1])); // was not from IOG
    int nsh = _nshl_vec[last_configuration - 1]; // was not from IOG
    if (i115 == 1) nsh += ies;
    str_name = "RCF_" + to_string(last_configuration); // was not from IOG
    str_type = "FLOAT64_(" + to_string(nsh) + ")";
    rec_name_list->append(str_name);
    rec_type_list->append(str_type);
    rec_ptr_list->append(&(_rcf[last_configuration - 1][0])); // was not from IOG
  }

  int dim3 = (idfc == 0) ? _number_of_scales : 1;
  int dc0m_size = 2 * dim3 * _number_of_spheres * _max_layers;
  double *dc0m = new double[dc0m_size]();
  int dc0_index = 0;
  for (int jxi468 = 1; jxi468 <= _number_of_scales; jxi468++) {
    if (_idfc != 0 && jxi468 > 1) continue;
    last_configuration = 0;
    for (int i162 = 1; i162 <= _number_of_spheres; i162++) {
      if (_iog_vec[i162 - 1] < i162) continue;
      else last_configuration++;
      int nsh = _nshl_vec[last_configuration - 1]; // was not from IOG
      int ici = (nsh + 1) / 2;
      if (i162 == 1) ici = ici + ies;
      for (int i157 = 0; i157 < ici; i157++) {
	double dc0_real, dc0_imag;
	dc0_real = real(_dc0_matrix[i157][i162 - 1][jxi468 - 1]);
	dc0_imag = imag(_dc0_matrix[i157][i162 - 1][jxi468 - 1]);
	dc0m[2 * dc0_index] = dc0_real;
	dc0m[2 * dc0_index + 1] = dc0_imag;
	dc0_index++;
      }
    }
  }
  str_type = "FLOAT64_(" + to_string(dc0m_size) + ")";
  rec_name_list->append("DC0M");
  rec_type_list->append(str_type);
  rec_ptr_list->append(dc0m);

  string *rec_names = rec_name_list->to_array();
  string *rec_types = rec_type_list->to_array();
  void **rec_pointers = rec_ptr_list->to_array();
  const int rec_num = rec_name_list->length();
  FileSchema *schema = new FileSchema(rec_num, rec_types, rec_names);
  HDFFile *hdf_file = HDFFile::from_schema(*schema, file_name, H5F_ACC_TRUNC);
  for (int ri = 0; ri < rec_num; ri++)
    hdf_file->write(rec_names[ri], rec_types[ri], rec_pointers[ri]);
  hdf_file->close();

  // Clean memory
  delete rec_name_list;
  delete rec_type_list;
  delete rec_ptr_list;
  delete[] dc0m;
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete schema;
  delete hdf_file;
}

void ScattererConfiguration::write_legacy(const std::string& file_name) {
  fstream output;
  int ies = (_use_external_sphere)? 1 : 0;
  int last_configuration;
  output.open(file_name.c_str(), ios::out | ios::binary);
  output.write(reinterpret_cast<char *>(&_number_of_spheres), sizeof(int));
  output.write(reinterpret_cast<char *>(&ies), sizeof(int));
  for (int i = 0; i < _number_of_spheres; i++)
    output.write(reinterpret_cast<char *>(&(_iog_vec[i])), sizeof(int));
  output.write(reinterpret_cast<char *>(&_exdc), sizeof(double));
  output.write(reinterpret_cast<char *>(&_wp), sizeof(double));
  output.write(reinterpret_cast<char *>(&_xip), sizeof(double));
  output.write(reinterpret_cast<char *>(&_idfc), sizeof(int));
  output.write(reinterpret_cast<char *>(&_number_of_scales), sizeof(int));
  for (int i = 0; i < _number_of_scales; i++)
    output.write(reinterpret_cast<char *>(&(_scale_vec[i])), sizeof(double));
  last_configuration = 0;
  for (int i115 = 1; i115 <= _number_of_spheres; i115++) {
    if (_iog_vec[i115 - 1] < i115) continue;
    else last_configuration++;
    output.write(reinterpret_cast<char *>(&(_nshl_vec[last_configuration - 1])), sizeof(int)); // was not from IOG
    output.write(reinterpret_cast<char *>(&(_radii_of_spheres[last_configuration - 1])), sizeof(double)); // was not from IOG
    int nsh = _nshl_vec[last_configuration - 1]; // was not from IOG
    if (i115 == 1) nsh += ies;
    for (int nsi = 0; nsi < nsh; nsi++)
      output.write(reinterpret_cast<char *>(&(_rcf[last_configuration - 1][nsi])), sizeof(double)); // was not from IOG
  }
  for (int jxi468 = 1; jxi468 <= _number_of_scales; jxi468++) {
    if (idfc != 0 && jxi468 > 1) continue;
    last_configuration = 0;
    for (int i162 = 1; i162 <= _number_of_spheres; i162++) {
      if (_iog_vec[i162 - 1] < i162) continue;
      else last_configuration++;
      int nsh = _nshl_vec[last_configuration - 1]; // was not from IOG
      int ici = (nsh + 1) / 2; // QUESTION: is integer division really intended here?
      if (i162 == 1) ici = ici + ies;
      for (int i157 = 0; i157 < ici; i157++) {
	double dc0_real, dc0_img;
	dc0_real = real(_dc0_matrix[i157][i162 - 1][jxi468 - 1]);
	dc0_img = imag(_dc0_matrix[i157][i162 - 1][jxi468 - 1]);
	// The FORTRAN code writes the complex numbers as a 16-byte long binary stream.
	// Here we assume that the 16 bytes are equally split in 8 bytes to represent the
	// real part and 8 bytes to represent the imaginary one.
	output.write(reinterpret_cast<char *>(&dc0_real), sizeof(double));
	output.write(reinterpret_cast<char *>(&dc0_img), sizeof(double));
      }
    }
  }
  output.close();
}

void ScattererConfiguration::write_formatted(const std::string& file_name) {
  const double evc = 6.5821188e-16;
  const double two_pi = acos(0.0) * 4.0;
  double *xi_vec;
  int ici;
  int last_configuration;
  int ies = (use_external_sphere)? 1: 0;
  FILE *output = fopen(file_name.c_str(), "w");
  int scale_type = -1;
  if (_reference_variable_name.compare("XIV") == 0) scale_type = 0;
  else if (_reference_variable_name.compare("WNS") == 0) scale_type = 1;
  else if (_reference_variable_name.compare("WLS") == 0) scale_type = 2;
  else if (_reference_variable_name.compare("PUS") == 0) scale_type = 3;
  else if (_reference_variable_name.compare("EVS") == 0) scale_type = 4;
  if (_idfc >= 0) { // Dielectric functions are constant or they depend on XI
    double  *pu_vec, *ev_vec, *wn_vec, *wl_vec;
    xi_vec = new double[_number_of_scales];
    pu_vec = new double[_number_of_scales];
    ev_vec = new double[_number_of_scales];
    wn_vec = new double[_number_of_scales];
    wl_vec = new double[_number_of_scales];
    switch (scale_type) {
    case 0:
      fprintf(output, "  JXI     XIV          WNS          WLS          PUS          EVS\n");
      for (int i = 0; i < _number_of_scales; i++) {
	xi_vec[i] = _scale_vec[i];
	pu_vec[i] = xi_vec[i] * wp;
	ev_vec[i] = pu_vec[i] * evc;
	wn_vec[i] = pu_vec[i] / 3.0e8;
	wl_vec[i] = two_pi / wn_vec[i];
	fprintf(
		output,
		"%5d%13.4lE%13.4lE%13.4lE%13.4lE%13.4lE\n",
		(i + 1),
		xi_vec[i],
		wn_vec[i],
		wl_vec[i],
		pu_vec[i],
		ev_vec[i]
		);
      }
      break;
    case 1:
      fprintf(output, "  JXI     WNS          WLS          PUS          EVS          XIV\n");
      for (int i = 0; i < _number_of_scales; i++) {
	xi_vec[i] = _scale_vec[i];
	pu_vec[i] = xi_vec[i] * wp;
	wn_vec[i] = pu_vec[i] / 3.0e8;
	wl_vec[i] = two_pi / wn_vec[i];
	ev_vec[i] = pu_vec[i] * evc;
	fprintf(
		output,
		"%5d%13.4lE%13.4lE%13.4lE%13.4lE%13.4lE\n",
		(i + 1),
		wn_vec[i],
		wl_vec[i],
		pu_vec[i],
		ev_vec[i],
		xi_vec[i]
		);
      }
      break;
    case 2:
      fprintf(output, "  JXI     WLS          WNS          PUS          EVS          XIV\n");
      for (int i = 0; i < _number_of_scales; i++) {
	xi_vec[i] = _scale_vec[i];
	pu_vec[i] = xi_vec[i] * wp;
	wn_vec[i] = pu_vec[i] / 3.0e8;
	wl_vec[i] = two_pi / wn_vec[i];
	ev_vec[i] = pu_vec[i] * evc;
	fprintf(
		output,
		"%5d%13.4lE%13.4lE%13.4lE%13.4lE%13.4lE\n",
		(i + 1),
		wl_vec[i],
		wn_vec[i],
		pu_vec[i],
		ev_vec[i],
		xi_vec[i]
		);
      }
      break;
    case 3:
      fprintf(output, "  JXI     PUS          WNS          WLS          EVS          XIV\n");
      for (int i = 0; i < _number_of_scales; i++) {
	xi_vec[i] = _scale_vec[i];
	pu_vec[i] = xi_vec[i] * wp;
	wn_vec[i] = pu_vec[i] / 3.0e8;
	wl_vec[i] = two_pi / wn_vec[i];
	ev_vec[i] = pu_vec[i] * evc;
	fprintf(
		output,
		"%5d%13.4lE%13.4lE%13.4lE%13.4lE%13.4lE\n",
		(i + 1),
		pu_vec[i],
		wn_vec[i],
		wl_vec[i],
		ev_vec[i],
		xi_vec[i]
		);
      }
      break;
    case 4:
      fprintf(output, "  JXI     EVS          WNS          WLS          PUS          XIV\n");
      for (int i = 0; i < _number_of_scales; i++) {
	xi_vec[i] = _scale_vec[i];
	pu_vec[i] = xi_vec[i] * wp;
	wn_vec[i] = pu_vec[i] / 3.0e8;
	wl_vec[i] = two_pi / wn_vec[i];
	ev_vec[i] = pu_vec[i] * evc;
	fprintf(
		output,
		"%5d%13.4lE%13.4lE%13.4lE%13.4lE%13.4lE\n",
		(i + 1),
		ev_vec[i],
		wn_vec[i],
		wl_vec[i],
		pu_vec[i],
		xi_vec[i]
		);
      }
      break;
    default:
      throw UnrecognizedConfigurationException(
					       "Wonrg parameter set: unrecognized scale variable type "
					       + reference_variable_name
					       );
      break;
    }
    // Clean memory
    delete[] xi_vec;
    delete[] pu_vec;
    delete[] ev_vec;
    delete[] wn_vec;
    delete[] wl_vec;
  } else { // idfc < 0, Dielectric functions are at XIP and XI is scale for dimensions
    double pu, wn;
    xi_vec = _scale_vec;
    pu = xip * wp;
    wn = pu / 3.0e8;
    fprintf(output, "          XIP          WN           WL           PU           EV\n");
    fprintf(output, "     %13.4lE", xip);
    fprintf(output, "%13.4lE", wn);
    fprintf(output, "%13.4lE", two_pi / wn);
    fprintf(output, "%13.4lE", pu);
    fprintf(output, "%13.4lE\n", pu * evc);
    fprintf(output, "  SCALE FACTORS XI\n");
    for (int i = 0; i < _number_of_scales; i++)
      fprintf(output, "%5d%13.4lE\n", (i + 1), xi_vec[i]);
  }
  if (_idfc != 0) {
    fprintf(output, "  DIELECTRIC CONSTANTS\n");
    last_configuration = 0;
    for (int i473 = 1; i473 <= _number_of_spheres; i473++) {
      if (_iog_vec[i473 - 1] != i473) continue;
      else last_configuration++;
      ici = (_nshl_vec[last_configuration - 1] + 1) / 2;
      if (i473 == 1) ici += ies;
      fprintf(output, " SPHERE N.%4d\n", i473);
      for (int ic472 = 0; ic472 < ici; ic472++) {
	double dc0_real = real(_dc0_matrix[ic472][i473 - 1][0]);
	double dc0_img = imag(_dc0_matrix[ic472][i473 - 1][0]);
	fprintf(output, "%5d %12.4lE%12.4lE\n", (ic472 + 1), dc0_real, dc0_img);
      }
    }
  } else {
    fprintf(output, "  DIELECTRIC FUNCTIONS\n");
    last_configuration = 0;
    for (int i478 = 1; i478 <= _number_of_spheres; i478++) {
      if (_iog_vec[i478 - 1] != i478) continue;
      last_configuration++;
      ici = (_nshl_vec[last_configuration - 1] + 1) / 2;
      if (i478 == 1) ici += ies;
      fprintf(output, " SPHERE N.%4d\n", i478);
      for (int ic477 = 1; ic477 <= ici; ic477++) {
	fprintf(output, " NONTRANSITION LAYER N.%2d, SCALE = %3s\n", ic477, _reference_variable_name.c_str());
	for (int jxi476 = 0; jxi476 < _number_of_scales; jxi476++) {
	  double dc0_real = real(_dc0_matrix[ic477 - 1][i478 - 1][jxi476]);
	  double dc0_img = imag(_dc0_matrix[ic477 - 1][i478 - 1][jxi476]);
	  fprintf(output, "%5d %12.4lE%12.4lE\n", (jxi476 + 1), dc0_real, dc0_img);
	}
      }
    }
  }
  fclose(output);
}

bool ScattererConfiguration::operator ==(const ScattererConfiguration &other) {
  if (_number_of_spheres != other._number_of_spheres) {
    return false;
  }
  if (_number_of_scales != other._number_of_scales) {
    return false;
  }
  if (_idfc != other._idfc) {
    return false;
  }
  if (_max_layers != other._max_layers) {
    return false;
  }
  if (_exdc != other._exdc) {
    return false;
  }
  if (_wp != other._wp) {
    return false;
  }
  if (_xip != other._xip) {
    return false;
  }
  //if (use_external_sphere != other.use_external_sphere) return false;
  for (int svi = 0; svi < _number_of_scales; svi++)
    if (_scale_vec[svi] != other._scale_vec[svi]) {
      return false;
    }
  if (_configurations != other._configurations) {
    return false;
  }
  int dim3 = (idfc == 0) ? _number_of_scales : 1;
  for (int ci = 1; ci <= _number_of_spheres; ci++) {
    if (_iog_vec[ci - 1] != other._iog_vec[ci - 1]) {
      return false;
    }
  }
  for (int ri = 0; ri < _configurations; ri++) {
    if (_radii_of_spheres[ri] != other._radii_of_spheres[ri]) {
      return false;
    }
    if (_nshl_vec[ri] != other._nshl_vec[ri]) {
      return false;
    }
    for (int rj = 0; rj < _nshl_vec[ri]; rj++) {
      if (_rcf[ri][rj] != other._rcf[ri][rj]) {
	return false;
      }
    } // rj loop
  } // ri loop
  for (int li = 0; li < _max_layers; li++) {
    for (int dj = 0; dj < _number_of_spheres; dj++) {
      for (int dk = 0; dk < dim3; dk++) {
	if (_dc0_matrix[li][dj][dk] != other._dc0_matrix[li][dj][dk]) {
	  return false;
	}
      } // dk loop
    } // dj loop
  }
  return true;
}
