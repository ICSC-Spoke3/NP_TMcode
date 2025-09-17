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

/*! \file outputs.cpp
 *
 * \brief Implementation of the code output format system.
 */
#include <cstdio>
#include <exception>
#include <string>
#include <cstring>
#include <hdf5.h>

#ifdef USE_MPI
#ifndef MPI_VERSION
#include <mpi.h>
#endif
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

using namespace std;

// >>> ClusterOutputInfo CLASS IMPLEMENTATION <<<
ClusterOutputInfo::ClusterOutputInfo(
  ScattererConfiguration *sc, GeometryConfiguration *gc,
  const mixMPI *mpidata, int first_xi, int xi_length
) {
  _skip_flag = 0;
  nsph = gc->number_of_spheres;
  li = gc->li;
  le = gc->le;
  lm = gc->l_max;
  mxndm = gc->mxndm;
  inpol = gc->in_pol;
  npnt = gc->npnt;
  npntts = gc->npntts;
  iavm = gc->iavm;
  isam = gc->isam;
  jwtm = gc->jwtm;
  // Get spherical constituent coordinates
  vec_x_coords = new double[nsph];
  vec_y_coords = new double[nsph];
  vec_z_coords = new double[nsph];
  for (int nsi = 0; nsi < nsph; nsi++) {
    vec_x_coords[nsi] = gc->get_sph_x(nsi);
    vec_y_coords[nsi] = gc->get_sph_y(nsi);
    vec_z_coords[nsi] = gc->get_sph_z(nsi);
  }
  // Get directional information
  th = gc->in_theta_start;
  thstp = gc->in_theta_step;
  thlst = gc->in_theta_end;
  ths = gc->sc_theta_start;
  thsstp = gc->sc_theta_step;
  thslst = gc->sc_theta_end;
  _num_theta = (thstp == 0.0) ? 1 : 1 + int((thlst - th) / thstp);
  _num_thetas = (thsstp == 0.0) ? 1 : 1 + int((thslst - ths) / thsstp);
  ph = gc->in_phi_start;
  phstp = gc->in_phi_step;
  phlst = gc->in_phi_end;
  phs = gc->sc_phi_start;
  phsstp = gc->sc_phi_step;
  phslst = gc->sc_phi_end;
  _num_phi = (phstp == 0.0) ? 1 : 1 + int((phlst - ph) / phstp);
  _num_phis = (phsstp == 0.0) ? 1 : 1 + int((phslst - phs) / phsstp);
  ndirs = _num_theta * _num_thetas * _num_phi * _num_phis;
  // Get scattering problem configuration
  double exdc = sc->exdc;
  exri = sqrt(exdc);
  idfc = sc->idfc;
  nxi = sc->number_of_scales;
  _first_xi = first_xi;
  xi_block_size = (xi_length == 0) ? nxi - (_first_xi - 1) : xi_length;
  vec_jxi = new int[xi_block_size]();
  vec_ier = new short[xi_block_size]();
  vec_vk = new double[xi_block_size]();
  vec_xi = new double[xi_block_size]();
  configurations = sc->configurations;
  vec_sphere_sizes = new double[xi_block_size * configurations]();
  vec_sphere_ref_indices = new dcomplex[xi_block_size * configurations]();
  vec_sphere_scs = new double[xi_block_size * configurations]();
  vec_sphere_abs = new double[xi_block_size * configurations]();
  vec_sphere_exs = new double[xi_block_size * configurations]();
  vec_sphere_albs = new double[xi_block_size * configurations]();
  vec_sphere_sqscs = new double[xi_block_size * configurations]();
  vec_sphere_sqabs = new double[xi_block_size * configurations]();
  vec_sphere_sqexs = new double[xi_block_size * configurations]();
  vec_fsas = new dcomplex[xi_block_size * configurations]();
  vec_qschus = new double[xi_block_size * configurations]();
  vec_pschus = new double[xi_block_size * configurations]();
  vec_s0mags = new double[xi_block_size * configurations]();
  vec_cosavs = new double[xi_block_size * configurations]();
  vec_raprs = new double[xi_block_size * configurations]();
  vec_tqek1 = new double[xi_block_size * configurations]();
  vec_tqsk1 = new double[xi_block_size * configurations]();
  vec_tqek2 = new double[xi_block_size * configurations]();
  vec_tqsk2 = new double[xi_block_size * configurations]();
  vec_fsat = new dcomplex[xi_block_size]();
  vec_qschut = new double[xi_block_size]();
  vec_pschut = new double[xi_block_size]();
  vec_s0magt = new double[xi_block_size]();
  tgs = 0.0;
  vec_scc1 = new double[xi_block_size]();
  vec_scc2 = new double[xi_block_size]();
  vec_abc1 = new double[xi_block_size]();
  vec_abc2 = new double[xi_block_size]();
  vec_exc1 = new double[xi_block_size]();
  vec_exc2 = new double[xi_block_size]();
  vec_albedc1 = new double[xi_block_size]();
  vec_albedc2 = new double[xi_block_size]();
  vec_qscamc1 = new double[xi_block_size]();
  vec_qscamc2 = new double[xi_block_size]();
  vec_qabsmc1 = new double[xi_block_size]();
  vec_qabsmc2 = new double[xi_block_size]();
  vec_qextmc1 = new double[xi_block_size]();
  vec_qextmc2 = new double[xi_block_size]();
  vec_sccrt1 = new double[xi_block_size]();
  vec_sccrt2 = new double[xi_block_size]();
  vec_abcrt1 = new double[xi_block_size]();
  vec_abcrt2 = new double[xi_block_size]();
  vec_excrt1 = new double[xi_block_size]();
  vec_excrt2 = new double[xi_block_size]();
  vec_fsac11 = new dcomplex[xi_block_size]();
  vec_fsac21 = new dcomplex[xi_block_size]();
  vec_fsac22 = new dcomplex[xi_block_size]();
  vec_fsac12 = new dcomplex[xi_block_size]();
  vec_qschuc1 = new double[xi_block_size]();
  vec_qschuc2 = new double[xi_block_size]();
  vec_pschuc1 = new double[xi_block_size]();
  vec_pschuc2 = new double[xi_block_size]();
  vec_s0magc1 = new double[xi_block_size]();
  vec_s0magc2 = new double[xi_block_size]();
  vec_cosavc1 = new double[xi_block_size]();
  vec_cosavc2 = new double[xi_block_size]();
  vec_raprc1 = new double[xi_block_size]();
  vec_raprc2 = new double[xi_block_size]();
  vec_fkc1 = new double[xi_block_size]();
  vec_fkc2 = new double[xi_block_size]();
  vec_dir_tidg = new double[_num_theta];
  vec_dir_pidg = new double[_num_phi];
  vec_dir_tsdg = new double[_num_thetas];
  vec_dir_psdg = new double[_num_phis];
  // Initialize directions (they are scale-independent)
  double cti = th, cpi = ph, cts = ths, cps = phs;
  for (int di = 0; di < _num_theta; di++) {
    vec_dir_tidg[di] = cti;
    cti += thstp;
  }
  for (int di = 0; di < _num_thetas; di++) {
    vec_dir_tsdg[di] = cts;
    cts += thsstp;
  }
  for (int di = 0; di < _num_phi; di++) {
    vec_dir_pidg[di] = cpi;
    cpi += phstp;
  }
  for (int di = 0; di < _num_phis; di++) {
    vec_dir_psdg[di] = cps;
    cps += phsstp;
  }
  vec_dir_scand = new double[ndirs]();
  vec_dir_cfmp = new double[ndirs]();
  vec_dir_sfmp = new double[ndirs]();
  vec_dir_cfsp = new double[ndirs]();
  vec_dir_sfsp = new double[ndirs]();
  vec_dir_un = new double[3 * ndirs]();
  vec_dir_uns = new double[3 * ndirs]();
  vec_dir_sas11 = new dcomplex[ndirs * configurations * xi_block_size]();
  vec_dir_sas21 = new dcomplex[ndirs * configurations * xi_block_size]();
  vec_dir_sas12 = new dcomplex[ndirs * configurations * xi_block_size]();
  vec_dir_sas22 = new dcomplex[ndirs * configurations * xi_block_size]();
  vec_dir_muls = new double[16 * ndirs * configurations * xi_block_size]();
  vec_dir_mulslr = new double[16 * ndirs * configurations * xi_block_size]();
  vec_dir_sat11 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sat21 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sat12 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sat22 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_scc1 = new double[ndirs * xi_block_size]();
  vec_dir_scc2 = new double[ndirs * xi_block_size]();
  vec_dir_abc1 = new double[ndirs * xi_block_size]();
  vec_dir_abc2 = new double[ndirs * xi_block_size]();
  vec_dir_exc1 = new double[ndirs * xi_block_size]();
  vec_dir_exc2 = new double[ndirs * xi_block_size]();
  vec_dir_albedc1 = new double[ndirs * xi_block_size]();
  vec_dir_albedc2 = new double[ndirs * xi_block_size]();
  vec_dir_qscc1 = new double[ndirs * xi_block_size]();
  vec_dir_qscc2 = new double[ndirs * xi_block_size]();
  vec_dir_qabc1 = new double[ndirs * xi_block_size]();
  vec_dir_qabc2 = new double[ndirs * xi_block_size]();
  vec_dir_qexc1 = new double[ndirs * xi_block_size]();
  vec_dir_qexc2 = new double[ndirs * xi_block_size]();
  vec_dir_sccrt1 = new double[ndirs * xi_block_size]();
  vec_dir_sccrt2 = new double[ndirs * xi_block_size]();
  vec_dir_abcrt1 = new double[ndirs * xi_block_size]();
  vec_dir_abcrt2 = new double[ndirs * xi_block_size]();
  vec_dir_excrt1 = new double[ndirs * xi_block_size]();
  vec_dir_excrt2 = new double[ndirs * xi_block_size]();
  vec_dir_fsac11 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_fsac21 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_fsac12 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_fsac22 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sac11 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sac21 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sac12 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sac22 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_qschuc1 = new double[ndirs * xi_block_size]();
  vec_dir_qschuc2 = new double[ndirs * xi_block_size]();
  vec_dir_pschuc1 = new double[ndirs * xi_block_size]();
  vec_dir_pschuc2 = new double[ndirs * xi_block_size]();
  vec_dir_s0magc1 = new double[ndirs * xi_block_size]();
  vec_dir_s0magc2 = new double[ndirs * xi_block_size]();
  vec_dir_cosavc1 = new double[ndirs * xi_block_size]();
  vec_dir_cosavc2 = new double[ndirs * xi_block_size]();
  vec_dir_raprc1 = new double[ndirs * xi_block_size]();
  vec_dir_raprc2 = new double[ndirs * xi_block_size]();
  vec_dir_flc1 = new double[ndirs * xi_block_size]();
  vec_dir_flc2 = new double[ndirs * xi_block_size]();
  vec_dir_frc1 = new double[ndirs * xi_block_size]();
  vec_dir_frc2 = new double[ndirs * xi_block_size]();
  vec_dir_fkc1 = new double[ndirs * xi_block_size]();
  vec_dir_fkc2 = new double[ndirs * xi_block_size]();
  vec_dir_fxc1 = new double[ndirs * xi_block_size]();
  vec_dir_fxc2 = new double[ndirs * xi_block_size]();
  vec_dir_fyc1 = new double[ndirs * xi_block_size]();
  vec_dir_fyc2 = new double[ndirs * xi_block_size]();
  vec_dir_fzc1 = new double[ndirs * xi_block_size]();
  vec_dir_fzc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqelc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqelc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqerc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqerc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqekc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqekc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqexc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqexc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqeyc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqeyc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqezc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqezc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqslc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqslc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsrc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsrc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqskc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqskc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsxc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsxc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsyc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsyc2 = new double[ndirs * xi_block_size]();
  vec_dir_tqszc1 = new double[ndirs * xi_block_size]();
  vec_dir_tqszc2 = new double[ndirs * xi_block_size]();
  vec_dir_mulc = new double[16 * ndirs * xi_block_size]();
  vec_dir_mulclr = new double[16 * ndirs * xi_block_size]();
}

ClusterOutputInfo::ClusterOutputInfo(const std::string &hdf5_name) {
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(hdf5_name, flags);
  herr_t status = hdf_file->get_status();
  string str_name, str_type;
  _skip_flag = 0;
  if (status == 0) {
    status = hdf_file->read("NSPH", "INT32_(1)", &nsph);
    status = hdf_file->read("LI", "INT32_(1)", &li);
    status = hdf_file->read("LE", "INT32_(1)", &le);
    status = hdf_file->read("LM", "INT32_(1)", &lm);
    long tmp;
    status = hdf_file->read("MXNDM", "INT64_(1)", &tmp);
    mxndm = (np_int)tmp;
    status = hdf_file->read("INPOL", "INT32_(1)", &inpol);
    status = hdf_file->read("NPNT", "INT32_(1)", &npnt);
    status = hdf_file->read("NPNTTS", "INT32_(1)", &npntts);
    status = hdf_file->read("IAVM", "INT32_(1)", &iavm);
    status = hdf_file->read("ISAM", "INT32_(1)", &isam);
    status = hdf_file->read("JWTM", "INT32_(1)", &jwtm);
    str_type = "FLOAT64_(" + to_string(nsph) + ")";
    vec_x_coords = new double[nsph];
    vec_y_coords = new double[nsph];
    vec_z_coords = new double[nsph];
    status = hdf_file->read("VEC_SPH_X", str_type, vec_x_coords);
    status = hdf_file->read("VEC_SPH_Y", str_type, vec_y_coords);
    status = hdf_file->read("VEC_SPH_Z", str_type, vec_z_coords);
    status = hdf_file->read("TH_START", "FLOAT64_(1)", &th);
    status = hdf_file->read("TH_STEP", "FLOAT64_(1)", &thstp);
    status = hdf_file->read("TH_END", "FLOAT64_(1)", &thlst);
    _num_theta = (thstp == 0.0) ? 1 : 1 + int((thlst - th) / thstp);
    status = hdf_file->read("THS_START", "FLOAT64_(1)", &ths);
    status = hdf_file->read("THS_STEP", "FLOAT64_(1)", &thsstp);
    status = hdf_file->read("THS_END", "FLOAT64_(1)", &thslst);
    _num_thetas = (thsstp == 0.0) ? 1 : 1 + int((thslst - ths) / thsstp);
    status = hdf_file->read("PH_START", "FLOAT64_(1)", &ph);
    status = hdf_file->read("PH_STEP", "FLOAT64_(1)", &phstp);
    status = hdf_file->read("PH_END", "FLOAT64_(1)", &phlst);
    _num_phi = (phstp == 0.0) ? 1 : 1 + int((phlst - ph) / phstp);
    status = hdf_file->read("PHS_START", "FLOAT64_(1)", &phs);
    status = hdf_file->read("PHS_STEP", "FLOAT64_(1)", &phsstp);
    status = hdf_file->read("PHS_END", "FLOAT64_(1)", &phslst);
    _num_phis = (phsstp == 0.0) ? 1 : 1 + int((phslst - phs) / phsstp);
    ndirs = _num_theta * _num_thetas * _num_phi * _num_phis;
    status = hdf_file->read("EXRI", "FLOAT64_(1)", &exri);
    status = hdf_file->read("IDFC", "INT32_(1)", &idfc);
    status = hdf_file->read("XI1", "INT32_(1)", &_first_xi);
    status = hdf_file->read("NXI", "INT32_(1)", &xi_block_size);
    nxi = (_first_xi == 1) ? xi_block_size : xi_block_size + _first_xi;
    str_type = "INT32_(" + to_string(xi_block_size) + ")";
    vec_jxi = new int[xi_block_size];
    status = hdf_file->read("VEC_JXI", str_type, vec_jxi);
    str_type = "INT16_(" + to_string(xi_block_size) + ")";
    vec_ier = new short[xi_block_size];
    status = hdf_file->read("VEC_IER", str_type, vec_ier);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    vec_vk = new double[xi_block_size];
    status = hdf_file->read("VEC_VK", str_type, vec_vk);
    vec_xi = new double[xi_block_size];
    status = hdf_file->read("VEC_VK", str_type, vec_xi);
    status = hdf_file->read("NCONF", "INT32_(1)", &configurations);
    str_type = "FLOAT64_(" + to_string(xi_block_size * configurations) + ")";
    vec_sphere_sizes = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_SIZES", str_type, vec_sphere_sizes);
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size * configurations) + ")";
    vec_sphere_ref_indices = new dcomplex[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_REFRI", str_type, vec_sphere_ref_indices);
    str_type = "FLOAT64_(" + to_string(xi_block_size * configurations) + ")";
    vec_sphere_scs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_SCS", str_type, vec_sphere_scs);
    vec_sphere_abs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_ABS", str_type, vec_sphere_abs);
    vec_sphere_exs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_EXS", str_type, vec_sphere_exs);
    vec_sphere_albs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_ALBS", str_type, vec_sphere_albs);
    vec_sphere_sqscs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_SQSCS", str_type, vec_sphere_sqscs);
    vec_sphere_sqabs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_SQABS", str_type, vec_sphere_sqabs);
    vec_sphere_sqexs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_SPH_SQEXS", str_type, vec_sphere_sqexs);
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size * configurations) + ")";
    vec_fsas = new dcomplex[xi_block_size * configurations];
    status = hdf_file->read("VEC_FSAS", str_type, vec_fsas);
    str_type = "FLOAT64_(" + to_string(xi_block_size * configurations) + ")";
    vec_qschus = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_QSCHUS", str_type, vec_qschus);
    vec_pschus = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_PSCHUS", str_type, vec_pschus);
    vec_s0mags = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_S0MAGS", str_type, vec_s0mags);
    vec_cosavs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_COSAVS", str_type, vec_cosavs);
    vec_raprs = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_RAPRS", str_type, vec_raprs);
    vec_tqek1 = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_TQEK1", str_type, vec_tqek1);
    vec_tqsk1 = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_TQSK1", str_type, vec_tqsk1);
    vec_tqek2 = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_TQEK2", str_type, vec_tqek2);
    vec_tqsk2 = new double[xi_block_size * configurations];
    status = hdf_file->read("VEC_TQSK2", str_type, vec_tqsk2);
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size) + ")";
    vec_fsat = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAT", str_type, vec_fsat);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    vec_qschut = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCHUT", str_type, vec_qschut);
    vec_pschut = new double[xi_block_size];
    status = hdf_file->read("VEC_PSCHUT", str_type, vec_pschut);
    vec_s0magt = new double[xi_block_size];
    status = hdf_file->read("VEC_S0MAGT", str_type, vec_s0magt);
    vec_scc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCC1", str_type, vec_scc1);
    vec_scc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCC2", str_type, vec_scc2);
    vec_abc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABC1", str_type, vec_abc1);
    vec_abc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABC2", str_type, vec_abc2);
    vec_exc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXC1", str_type, vec_exc1);
    vec_exc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXC2", str_type, vec_exc2);
    vec_albedc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_ALBEDC1", str_type, vec_albedc1);
    vec_albedc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_ALBEDC2", str_type, vec_albedc2);
    vec_qscamc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCAMC1", str_type, vec_qscamc1);
    vec_qscamc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCAMC2", str_type, vec_qscamc2);
    vec_qabsmc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_QABSMC1", str_type, vec_qabsmc1);
    vec_qabsmc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_QABSMC2", str_type, vec_qabsmc2);
    vec_qextmc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_QEXTMC1", str_type, vec_qextmc1);
    vec_qextmc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_QEXTMC2", str_type, vec_qextmc2);
    vec_sccrt1 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCCRT1", str_type, vec_sccrt1);
    vec_sccrt2 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCCRT2", str_type, vec_sccrt2);
    vec_abcrt1 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABCRT1", str_type, vec_abcrt1);
    vec_abcrt2 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABCRT2", str_type, vec_abcrt2);
    vec_excrt1 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXCRT1", str_type, vec_excrt1);
    vec_excrt2 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXCRT2", str_type, vec_excrt2);
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size) + ")";
    vec_fsac11 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAC11", str_type, vec_fsac11);
    vec_fsac21 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAC21", str_type, vec_fsac21);
    vec_fsac22 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAC22", str_type, vec_fsac22);
    vec_fsac12 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAC12", str_type, vec_fsac12);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    vec_qschuc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCHUC1", str_type, vec_qschuc1);
    vec_qschuc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCHUC2", str_type, vec_qschuc2);
    vec_pschuc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_PSCHUC1", str_type, vec_pschuc1);
    vec_pschuc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_PSCHUC2", str_type, vec_pschuc2);
    vec_s0magc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_S0MAGC1", str_type, vec_s0magc1);
    vec_s0magc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_S0MAGC2", str_type, vec_s0magc2);
    vec_cosavc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_COSAVC1", str_type, vec_cosavc1);
    vec_cosavc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_COSAVC2", str_type, vec_cosavc2);
    vec_raprc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_RAPRC1", str_type, vec_raprc1);
    vec_raprc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_RAPRC2", str_type, vec_raprc2);
    vec_fkc1 = new double[xi_block_size];
    status = hdf_file->read("VEC_FKC1", str_type, vec_fkc1);
    vec_fkc2 = new double[xi_block_size];
    status = hdf_file->read("VEC_FKC2", str_type, vec_fkc2);
    // Initialize directions (they are scale-independent)
    vec_dir_tidg = new double[_num_theta];
    vec_dir_pidg = new double[_num_phi];
    vec_dir_tsdg = new double[_num_thetas];
    vec_dir_psdg = new double[_num_phis];
    double cti = th, cpi = ph, cts = ths, cps = phs;
    for (int di = 0; di < _num_theta; di++) {
      vec_dir_tidg[di] = cti;
      cti += thstp;
    }
    for (int di = 0; di < _num_thetas; di++) {
      vec_dir_tsdg[di] = cts;
      cts += thsstp;
    }
    for (int di = 0; di < _num_phi; di++) {
      vec_dir_pidg[di] = cpi;
      cpi += phstp;
    }
    for (int di = 0; di < _num_phis; di++) {
      vec_dir_psdg[di] = cps;
      cps += phsstp;
    }
    str_type = "FLOAT64_(" + to_string(ndirs) + ")";
    vec_dir_scand = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SCAN", str_type, vec_dir_scand);
    vec_dir_cfmp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_CFMP", str_type, vec_dir_cfmp);
    vec_dir_cfsp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_CFSP", str_type, vec_dir_cfsp);
    vec_dir_sfmp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SFMP", str_type, vec_dir_sfmp);
    vec_dir_sfsp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SFSP", str_type, vec_dir_sfsp);
    str_type = "FLOAT64_(" + to_string(3 * ndirs) + ")";
    vec_dir_un = new double[3 * ndirs];
    status = hdf_file->read("VEC_DIR_UN", str_type, vec_dir_un);
    vec_dir_uns = new double[3 * ndirs];
    status = hdf_file->read("VEC_DIR_UNS", str_type, vec_dir_uns);
    str_type = "FLOAT64_(" + to_string(2 * ndirs * configurations * xi_block_size) + ")";
    vec_dir_sas11 = new dcomplex[ndirs * configurations * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS11", str_type, vec_dir_sas11);
    vec_dir_sas21 = new dcomplex[ndirs * configurations * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS21", str_type, vec_dir_sas21);
    vec_dir_sas12 = new dcomplex[ndirs * configurations * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS12", str_type, vec_dir_sas12);
    vec_dir_sas22 = new dcomplex[ndirs * configurations * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS22", str_type, vec_dir_sas22);
    str_type = "FLOAT64_(" + to_string(16 * ndirs * configurations * xi_block_size) + ")";
    vec_dir_muls = new double[16 * ndirs *configurations * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULS", str_type, vec_dir_muls);
    vec_dir_mulslr = new double[16 * ndirs *configurations * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULSLR", str_type, vec_dir_mulslr);
    str_type = "FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")";
    vec_dir_sat11 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAT11", str_type, vec_dir_sat11);
    vec_dir_sat21 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAT21", str_type, vec_dir_sat21);
    vec_dir_sat12 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAT12", str_type, vec_dir_sat12);
    vec_dir_sat22 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAT22", str_type, vec_dir_sat22);
    str_type = "FLOAT64_(" + to_string(ndirs * xi_block_size) + ")";
    vec_dir_scc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCC1", str_type, vec_dir_scc1);
    vec_dir_scc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCC2", str_type, vec_dir_scc2);
    vec_dir_abc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABC1", str_type, vec_dir_abc1);
    vec_dir_abc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABC2", str_type, vec_dir_abc2);
    vec_dir_exc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXC1", str_type, vec_dir_exc1);
    vec_dir_exc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXC2", str_type, vec_dir_exc2);
    vec_dir_albedc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ALBEDC1", str_type, vec_dir_albedc1);
    vec_dir_albedc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ALBEDC2", str_type, vec_dir_albedc2);
    vec_dir_qscc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QSCC1", str_type, vec_dir_qscc1);
    vec_dir_qscc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QSCC2", str_type, vec_dir_qscc2);
    vec_dir_qabc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QABC1", str_type, vec_dir_qabc1);
    vec_dir_qabc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QABC2", str_type, vec_dir_qabc2);
    vec_dir_qexc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QEXC1", str_type, vec_dir_qexc1);
    vec_dir_qexc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QEXC2", str_type, vec_dir_qexc2);
    vec_dir_sccrt1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCCRT1", str_type, vec_dir_sccrt1);
    vec_dir_sccrt2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCCRT2", str_type, vec_dir_sccrt2);
    vec_dir_abcrt1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABCRT1", str_type, vec_dir_abcrt1);
    vec_dir_abcrt2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABCRT2", str_type, vec_dir_abcrt2);
    vec_dir_excrt1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXCRT1", str_type, vec_dir_excrt1);
    vec_dir_excrt2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXCRT2", str_type, vec_dir_excrt2);
    str_type = "FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")";
    vec_dir_fsac11 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAC11", str_type, vec_dir_fsac11);
    vec_dir_fsac21 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAC21", str_type, vec_dir_fsac21);
    vec_dir_fsac12 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAC12", str_type, vec_dir_fsac12);
    vec_dir_fsac22 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAC22", str_type, vec_dir_fsac22);
    vec_dir_sac11 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAC11", str_type, vec_dir_sac11);
    vec_dir_sac21 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAC21", str_type, vec_dir_sac21);
    vec_dir_sac12 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAC12", str_type, vec_dir_sac12);
    vec_dir_sac22 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAC22", str_type, vec_dir_sac22);
    str_type = "FLOAT64_(" + to_string(ndirs * xi_block_size) + ")";
    vec_dir_qschuc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QSCHUC1", str_type, vec_dir_qschuc1);
    vec_dir_qschuc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QSCHUC2", str_type, vec_dir_qschuc2);
    vec_dir_pschuc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_PSCHUC1", str_type, vec_dir_pschuc1);
    vec_dir_pschuc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_PSCHUC2", str_type, vec_dir_pschuc2);
    vec_dir_s0magc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_S0MAGC1", str_type, vec_dir_s0magc1);
    vec_dir_s0magc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_S0MAGC2", str_type, vec_dir_s0magc2);
    vec_dir_cosavc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_COSAVC1", str_type, vec_dir_cosavc1);
    vec_dir_cosavc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_COSAVC2", str_type, vec_dir_cosavc2);
    vec_dir_raprc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_RAPRC1", str_type, vec_dir_raprc1);
    vec_dir_raprc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_RAPRC2", str_type, vec_dir_raprc2);
    vec_dir_flc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FLC1", str_type, vec_dir_flc1);
    vec_dir_flc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FLC2", str_type, vec_dir_flc2);
    vec_dir_frc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FRC1", str_type, vec_dir_frc1);
    vec_dir_frc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FRC2", str_type, vec_dir_frc2);
    vec_dir_fkc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FKC1", str_type, vec_dir_fkc1);
    vec_dir_fkc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FKC2", str_type, vec_dir_fkc2);
    vec_dir_fxc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FXC1", str_type, vec_dir_fxc1);
    vec_dir_fxc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FXC2", str_type, vec_dir_fxc2);
    vec_dir_fyc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FYC1", str_type, vec_dir_fyc1);
    vec_dir_fyc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FYC2", str_type, vec_dir_fyc2);
    vec_dir_fzc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FZC1", str_type, vec_dir_fzc1);
    vec_dir_fzc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FZC2", str_type, vec_dir_fzc2);
    vec_dir_tqelc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQELC1", str_type, vec_dir_tqelc1);
    vec_dir_tqelc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQELC2", str_type, vec_dir_tqelc2);
    vec_dir_tqerc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQERC1", str_type, vec_dir_tqerc1);
    vec_dir_tqerc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQERC2", str_type, vec_dir_tqerc2);
    vec_dir_tqekc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEKC1", str_type, vec_dir_tqekc1);
    vec_dir_tqekc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEKC2", str_type, vec_dir_tqekc2);
    vec_dir_tqexc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEXC1", str_type, vec_dir_tqexc1);
    vec_dir_tqexc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEXC2", str_type, vec_dir_tqexc2);
    vec_dir_tqeyc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEYC1", str_type, vec_dir_tqeyc1);
    vec_dir_tqeyc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEYC2", str_type, vec_dir_tqeyc2);
    vec_dir_tqezc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEZC1", str_type, vec_dir_tqezc1);
    vec_dir_tqezc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEZC2", str_type, vec_dir_tqezc2);
    vec_dir_tqslc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSLC1", str_type, vec_dir_tqslc1);
    vec_dir_tqslc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSLC2", str_type, vec_dir_tqslc2);
    vec_dir_tqsrc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSRC1", str_type, vec_dir_tqsrc1);
    vec_dir_tqsrc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSRC2", str_type, vec_dir_tqsrc2);
    vec_dir_tqskc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSKC1", str_type, vec_dir_tqskc1);
    vec_dir_tqskc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSKC2", str_type, vec_dir_tqskc2);
    vec_dir_tqsxc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSXC1", str_type, vec_dir_tqsxc1);
    vec_dir_tqsxc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSXC2", str_type, vec_dir_tqsxc2);
    vec_dir_tqsyc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSYC1", str_type, vec_dir_tqsyc1);
    vec_dir_tqsyc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSYC2", str_type, vec_dir_tqsyc2);
    vec_dir_tqszc1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSZC1", str_type, vec_dir_tqszc1);
    vec_dir_tqszc2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSZC2", str_type, vec_dir_tqszc2);
    str_type = "FLOAT64_(" + to_string(16 * ndirs * xi_block_size) + ")";
    vec_dir_mulc = new double[16 * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULC", str_type, vec_dir_mulc);
    vec_dir_mulclr = new double[16 * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULCLR", str_type, vec_dir_mulclr);
    status = hdf_file->close();
    delete hdf_file;
  } else {
    if (hdf_file != NULL) delete hdf_file;
    UnrecognizedFormatException ex("Error: " + hdf5_name + " not recognized as a valid HDF5 file!");
    throw ex;
  }
}

ClusterOutputInfo::ClusterOutputInfo(const int flag) {
  /*
    create a dummy placeholder just to know I should skip MPI_Send and MPI_Recv
  */
  if (flag == 1) {
    _skip_flag = 1;
  } else {
    UnrecognizedOutputInfo ex(flag);
    throw ex;
  }
}

ClusterOutputInfo::~ClusterOutputInfo() {
  if (_skip_flag != 1) {
    delete[] vec_x_coords;
    delete[] vec_y_coords;
    delete[] vec_z_coords;
    delete[] vec_jxi;
    delete[] vec_ier;
    delete[] vec_vk;
    delete[] vec_xi;
    delete[] vec_sphere_sizes;
    delete[] vec_sphere_ref_indices;
    delete[] vec_sphere_scs;
    delete[] vec_sphere_abs;
    delete[] vec_sphere_exs;
    delete[] vec_sphere_albs;
    delete[] vec_sphere_sqscs;
    delete[] vec_sphere_sqabs;
    delete[] vec_sphere_sqexs;
    delete[] vec_fsas;
    delete[] vec_qschus;
    delete[] vec_pschus;
    delete[] vec_s0mags;
    delete[] vec_cosavs;
    delete[] vec_raprs;
    delete[] vec_tqek1;
    delete[] vec_tqsk1;
    delete[] vec_tqek2;
    delete[] vec_tqsk2;
    delete[] vec_fsat;
    delete[] vec_qschut;
    delete[] vec_pschut;
    delete[] vec_s0magt;
    delete[] vec_scc1;
    delete[] vec_scc2;
    delete[] vec_abc1;
    delete[] vec_abc2;
    delete[] vec_exc1;
    delete[] vec_exc2;
    delete[] vec_albedc1;
    delete[] vec_albedc2;
    delete[] vec_sccrt1;
    delete[] vec_sccrt2;
    delete[] vec_abcrt1;
    delete[] vec_abcrt2;
    delete[] vec_excrt1;
    delete[] vec_excrt2;
    delete[] vec_fsac11;
    delete[] vec_fsac21;
    delete[] vec_fsac22;
    delete[] vec_fsac12;
    delete[] vec_qschuc1;
    delete[] vec_qschuc2;
    delete[] vec_pschuc1;
    delete[] vec_pschuc2;
    delete[] vec_s0magc1;
    delete[] vec_s0magc2;
    delete[] vec_cosavc1;
    delete[] vec_cosavc2;
    delete[] vec_raprc1;
    delete[] vec_raprc2;
    delete[] vec_fkc1;
    delete[] vec_fkc2;
    delete[] vec_dir_tidg;
    delete[] vec_dir_pidg;
    delete[] vec_dir_tsdg;
    delete[] vec_dir_psdg;
    delete[] vec_dir_scand;
    delete[] vec_dir_cfmp;
    delete[] vec_dir_sfmp;
    delete[] vec_dir_cfsp;
    delete[] vec_dir_sfsp;
    delete[] vec_dir_un;
    delete[] vec_dir_uns;
    delete[] vec_dir_sas11;
    delete[] vec_dir_sas21;
    delete[] vec_dir_sas12;
    delete[] vec_dir_sas22;
    delete[] vec_dir_muls;
    delete[] vec_dir_mulslr;
    delete[] vec_dir_sat11;
    delete[] vec_dir_sat21;
    delete[] vec_dir_sat12;
    delete[] vec_dir_sat22;
    delete[] vec_dir_scc1;
    delete[] vec_dir_scc2;
    delete[] vec_dir_abc1;
    delete[] vec_dir_abc2;
    delete[] vec_dir_exc1;
    delete[] vec_dir_exc2;
    delete[] vec_dir_albedc1;
    delete[] vec_dir_albedc2;
    delete[] vec_dir_qscc1;
    delete[] vec_dir_qscc2;
    delete[] vec_dir_qabc1;
    delete[] vec_dir_qabc2;
    delete[] vec_dir_qexc1;
    delete[] vec_dir_qexc2;
    delete[] vec_qscamc1;
    delete[] vec_qscamc2;
    delete[] vec_qabsmc1;
    delete[] vec_qabsmc2;
    delete[] vec_qextmc1;
    delete[] vec_qextmc2;
    delete[] vec_dir_sccrt1;
    delete[] vec_dir_sccrt2;
    delete[] vec_dir_abcrt1;
    delete[] vec_dir_abcrt2;
    delete[] vec_dir_excrt1;
    delete[] vec_dir_excrt2;
    delete[] vec_dir_fsac11;
    delete[] vec_dir_fsac21;
    delete[] vec_dir_fsac12;
    delete[] vec_dir_fsac22;
    delete[] vec_dir_sac11;
    delete[] vec_dir_sac21;
    delete[] vec_dir_sac12;
    delete[] vec_dir_sac22;
    delete[] vec_dir_qschuc1;
    delete[] vec_dir_qschuc2;
    delete[] vec_dir_pschuc1;
    delete[] vec_dir_pschuc2;
    delete[] vec_dir_s0magc1;
    delete[] vec_dir_s0magc2;
    delete[] vec_dir_cosavc1;
    delete[] vec_dir_cosavc2;
    delete[] vec_dir_raprc1;
    delete[] vec_dir_raprc2;
    delete[] vec_dir_flc1;
    delete[] vec_dir_flc2;
    delete[] vec_dir_frc1;
    delete[] vec_dir_frc2;
    delete[] vec_dir_fkc1;
    delete[] vec_dir_fkc2;
    delete[] vec_dir_fxc1;
    delete[] vec_dir_fxc2;
    delete[] vec_dir_fyc1;
    delete[] vec_dir_fyc2;
    delete[] vec_dir_fzc1;
    delete[] vec_dir_fzc2;
    delete[] vec_dir_tqelc1;
    delete[] vec_dir_tqelc2;
    delete[] vec_dir_tqerc1;
    delete[] vec_dir_tqerc2;
    delete[] vec_dir_tqekc1;
    delete[] vec_dir_tqekc2;
    delete[] vec_dir_tqexc1;
    delete[] vec_dir_tqexc2;
    delete[] vec_dir_tqeyc1;
    delete[] vec_dir_tqeyc2;
    delete[] vec_dir_tqezc1;
    delete[] vec_dir_tqezc2;
    delete[] vec_dir_tqslc1;
    delete[] vec_dir_tqslc2;
    delete[] vec_dir_tqsrc1;
    delete[] vec_dir_tqsrc2;
    delete[] vec_dir_tqskc1;
    delete[] vec_dir_tqskc2;
    delete[] vec_dir_tqsxc1;
    delete[] vec_dir_tqsxc2;
    delete[] vec_dir_tqsyc1;
    delete[] vec_dir_tqsyc2;
    delete[] vec_dir_tqszc1;
    delete[] vec_dir_tqszc2;
    delete[] vec_dir_mulc;
    delete[] vec_dir_mulclr;
  }
}

long ClusterOutputInfo::compute_size(
  ScattererConfiguration *sc, GeometryConfiguration *gc,
  int first_xi, int xi_length
) {
  long result = sizeof(np_int);
  result += 21 * sizeof(int);
  result += 14 * sizeof(double);
  result += 121 * sizeof(long);
  int _nsph = gc->number_of_spheres;
  double _th = gc->in_theta_start;
  double _thstp = gc->in_theta_step;
  double _thlst = gc->in_theta_end;
  double _ths = gc->sc_theta_start;
  double _thsstp = gc->sc_theta_step;
  double _thslst = gc->sc_theta_end;
  int num_theta = 1 + int((_thlst - _th) / _thstp);
  int num_thetas = 1 + int((_thslst - _ths) / _thsstp);
  double _ph = gc->in_phi_start;
  double _phstp = gc->in_phi_step;
  double _phlst = gc->in_phi_end;
  double _phs = gc->sc_phi_start;
  double _phsstp = gc->sc_phi_step;
  double _phslst = gc->sc_phi_end;
  int num_phi = 1 + int((_phlst - _ph) / _phstp);
  int num_phis = 1 + int((_phslst - _phs) / _phsstp);
  int _ndirs = num_theta * num_thetas * num_phi * num_phis;
  int _nxi = sc->number_of_scales;
  int _xi_block_size = (xi_length == 0) ? _nxi : xi_length;
  int _configurations = sc->configurations;
  result += 3 * _nsph * sizeof(double); // sphere coordinate vectors
  result += _xi_block_size * sizeof(int); // scale index vector
  result += _xi_block_size * sizeof(short); // error code vector
  result += 2 * _xi_block_size * sizeof(double); // scale vectors
  result += _xi_block_size * _configurations * sizeof(dcomplex); // refraction indices vector
  result += 5 * _xi_block_size * _configurations * sizeof(double); // sphere sizes, albedos and cross-sections
  result += 3 * _xi_block_size * _configurations * sizeof(double); // cross-section to geometric section ratios
  result += _xi_block_size * _configurations * sizeof(dcomplex); // fsas vector
  result += 9 * _xi_block_size * _configurations * sizeof(double); // up to tqsk2 vector
  result += _xi_block_size * sizeof(dcomplex); // fsat vector
  result += 3 * _xi_block_size * sizeof(double); // up to s0magt vector
  result += 20 * _xi_block_size * sizeof(double); // up to excrtt vector
  result += 4 * _xi_block_size * sizeof(dcomplex); // up to fsac12 vector
  result += 12 * _xi_block_size * sizeof(double); // up to fkc vector
  result += num_theta * sizeof(double); // vec_dir_tidg;
  result += num_thetas * sizeof(double); // vec_dir_tsdg;
  result += num_phi * sizeof(double); // vec_dir_pidg;
  result += num_phis * sizeof(double); // vec_dir_psdg;
  result += 11 * _ndirs * sizeof(double); // up to dir_uns vector
  result += 4 * _ndirs * _configurations * _xi_block_size * sizeof(dcomplex); // up to dir_sas22 vector
  result += 32 * _ndirs * _configurations * _xi_block_size * sizeof(double); // up to dir_mulslr vector
  result += 4 * _ndirs * _xi_block_size * sizeof(dcomplex); // up to dir_sat22 vector
  result += 20 * _ndirs * _xi_block_size * sizeof(double); // up to dir_excrt vector
  result += 8 * _ndirs * _xi_block_size * sizeof(dcomplex); // up to dir_sac22 vector
  result += 80 * _ndirs * _xi_block_size * sizeof(double); // up to dir_mulclr vector
  
  return result;
}

long ClusterOutputInfo::compute_size() {
  long result = sizeof(np_int);
  result += 21 * sizeof(int);
  result += 14 * sizeof(double);
  result += 121 * sizeof(long);
  result += 3 * nsph * sizeof(double); // sphere coordinate vectors
  result += xi_block_size * sizeof(int); // scale index vector
  result += xi_block_size * sizeof(short); // error code vector
  result += 2 * xi_block_size * sizeof(double); // scale vectors
  result += xi_block_size * configurations * sizeof(dcomplex); // refraction indices vector
  result += 5 * xi_block_size * configurations * sizeof(double); // sphere sizes, albedos and cross-sections
  result += 3 * xi_block_size * configurations * sizeof(double); // cross-sections to geometric section ratios
  result += configurations * sizeof(double); // sphere geometric sections
  result += xi_block_size * configurations * sizeof(dcomplex); // fsas vector
  result += 9 * xi_block_size * configurations * sizeof(double); // up to tqsk2 vector
  result += xi_block_size * sizeof(dcomplex); // fsat vector
  result += 3 * xi_block_size * sizeof(double); // up to s0magt vector
  result += 20 * xi_block_size * sizeof(double); // up to excrtt vector
  result += 4 * xi_block_size * sizeof(dcomplex); // up to fsac12 vector
  result += 12 * xi_block_size * sizeof(double); // up to fkc vector
  result += _num_theta * sizeof(double); // vec_dir_tidg;
  result += _num_thetas * sizeof(double); // vec_dir_tsdg;
  result += _num_phi * sizeof(double); // vec_dir_pidg;
  result += _num_phis * sizeof(double); // vec_dir_psdg;
  result += 11 * ndirs * sizeof(double); // up to dir_uns vector
  result += 4 * ndirs * configurations * xi_block_size * sizeof(dcomplex); // up to dir_sas22 vector
  result += 32 * ndirs * configurations * xi_block_size * sizeof(double); // up to dir_mulslr vector
  result += 4 * ndirs * xi_block_size * sizeof(dcomplex); // up to dir_sat22 vector
  result += 20 * ndirs * xi_block_size * sizeof(double); // up to dir_excrt vector
  result += 8 * ndirs * xi_block_size * sizeof(dcomplex); // up to dir_sac22 vector
  result += 80 * ndirs * xi_block_size * sizeof(double); // up to dir_mulclr vector
  return result;
}

int ClusterOutputInfo::insert(const ClusterOutputInfo &rhs) {
  int result = 0;
  if (rhs.skip_flag != 1) {
    result += (rhs.nsph == nsph) ? 0 : 1;
    result += (rhs.inpol == inpol) ? 0 : 1;
    result += (rhs.iavm == iavm) ? 0 : 1;
    result += (rhs.isam == isam) ? 0 : 1;
    result += (rhs._num_theta == _num_theta) ? 0 : 1;
    result += (rhs._num_thetas == _num_thetas) ? 0 : 1;
    result += (rhs._num_phi == _num_phi) ? 0 : 1;
    result += (rhs._num_phis == _num_phis) ? 0 : 1;
    result += (rhs.ndirs == ndirs) ? 0 : 1;
    result += (rhs.exri == exri) ? 0 : 1;
    result += (rhs.idfc == idfc) ? 0 : 1;
    result += (rhs.configurations == configurations) ? 0 : 1;
    if (result == 0) {
      int offset, chunk_size, xi1;
      xi1 = rhs._first_xi;
      // Insert vectors whose size depends on wavelengths
      offset = xi1 - _first_xi;
      chunk_size = rhs.xi_block_size;
      memcpy(vec_jxi + offset, rhs.vec_jxi, chunk_size * sizeof(int));
      memcpy(vec_ier + offset, rhs.vec_ier, chunk_size * sizeof(short));
      memcpy(vec_vk + offset, rhs.vec_vk, chunk_size * sizeof(double));
      memcpy(vec_xi + offset, rhs.vec_xi, chunk_size * sizeof(double));
      memcpy(vec_fsat + offset, rhs.vec_fsat, chunk_size * sizeof(dcomplex));
      memcpy(vec_qschut + offset, rhs.vec_qschut, chunk_size * sizeof(double));
      memcpy(vec_pschut + offset, rhs.vec_pschut, chunk_size * sizeof(double));
      memcpy(vec_s0magt + offset, rhs.vec_s0magt, chunk_size * sizeof(double));
      memcpy(vec_scc1 + offset, rhs.vec_scc1, chunk_size * sizeof(double));
      memcpy(vec_scc2 + offset, rhs.vec_scc2, chunk_size * sizeof(double));
      memcpy(vec_abc1 + offset, rhs.vec_abc1, chunk_size * sizeof(double));
      memcpy(vec_abc2 + offset, rhs.vec_abc2, chunk_size * sizeof(double));
      memcpy(vec_exc1 + offset, rhs.vec_exc1, chunk_size * sizeof(double));
      memcpy(vec_exc2 + offset, rhs.vec_exc2, chunk_size * sizeof(double));
      memcpy(vec_albedc1 + offset, rhs.vec_albedc1, chunk_size * sizeof(double));
      memcpy(vec_albedc2 + offset, rhs.vec_albedc2, chunk_size * sizeof(double));
      memcpy(vec_qscamc1 + offset, rhs.vec_qscamc1, chunk_size * sizeof(double));
      memcpy(vec_qscamc2 + offset, rhs.vec_qscamc2, chunk_size * sizeof(double));
      memcpy(vec_qabsmc1 + offset, rhs.vec_qabsmc1, chunk_size * sizeof(double));
      memcpy(vec_qabsmc2 + offset, rhs.vec_qabsmc2, chunk_size * sizeof(double));
      memcpy(vec_qextmc1 + offset, rhs.vec_qextmc1, chunk_size * sizeof(double));
      memcpy(vec_qextmc2 + offset, rhs.vec_qextmc2, chunk_size * sizeof(double));
      memcpy(vec_sccrt1 + offset, rhs.vec_sccrt1, chunk_size * sizeof(double));
      memcpy(vec_sccrt2 + offset, rhs.vec_sccrt2, chunk_size * sizeof(double));
      memcpy(vec_abcrt1 + offset, rhs.vec_abcrt1, chunk_size * sizeof(double));
      memcpy(vec_abcrt2 + offset, rhs.vec_abcrt2, chunk_size * sizeof(double));
      memcpy(vec_excrt1 + offset, rhs.vec_excrt1, chunk_size * sizeof(double));
      memcpy(vec_excrt2 + offset, rhs.vec_excrt2, chunk_size * sizeof(double));
      memcpy(vec_fsac11 + offset, rhs.vec_fsac11, chunk_size * sizeof(dcomplex));
      memcpy(vec_fsac21 + offset, rhs.vec_fsac21, chunk_size * sizeof(dcomplex));
      memcpy(vec_fsac22 + offset, rhs.vec_fsac22, chunk_size * sizeof(dcomplex));
      memcpy(vec_fsac12 + offset, rhs.vec_fsac12, chunk_size * sizeof(dcomplex));
      memcpy(vec_qschuc1 + offset, rhs.vec_qschuc1, chunk_size * sizeof(double));
      memcpy(vec_qschuc2 + offset, rhs.vec_qschuc2, chunk_size * sizeof(double));
      memcpy(vec_pschuc1 + offset, rhs.vec_pschuc1, chunk_size * sizeof(double));
      memcpy(vec_pschuc2 + offset, rhs.vec_pschuc2, chunk_size * sizeof(double));
      memcpy(vec_s0magc1 + offset, rhs.vec_s0magc1, chunk_size * sizeof(double));
      memcpy(vec_s0magc2 + offset, rhs.vec_s0magc2, chunk_size * sizeof(double));
      memcpy(vec_cosavc1 + offset, rhs.vec_cosavc1, chunk_size * sizeof(double));
      memcpy(vec_cosavc2 + offset, rhs.vec_cosavc2, chunk_size * sizeof(double));
      memcpy(vec_raprc1 + offset, rhs.vec_raprc1, chunk_size * sizeof(double));
      memcpy(vec_raprc2 + offset, rhs.vec_raprc2, chunk_size * sizeof(double));
      memcpy(vec_fkc1 + offset, rhs.vec_fkc1, chunk_size * sizeof(double));
      memcpy(vec_fkc2 + offset, rhs.vec_fkc2, chunk_size * sizeof(double));
      // Insert vectors of multiple configuration values per scale
      offset = (xi1 - _first_xi) * configurations;
      chunk_size = rhs.xi_block_size * configurations;
      memcpy(vec_sphere_sizes + offset, rhs.vec_sphere_sizes, chunk_size * sizeof(double));
      memcpy(vec_sphere_ref_indices + offset, rhs.vec_sphere_ref_indices, chunk_size * sizeof(dcomplex));
      memcpy(vec_sphere_scs + offset, rhs.vec_sphere_scs, chunk_size * sizeof(double));
      memcpy(vec_sphere_abs + offset, rhs.vec_sphere_abs, chunk_size * sizeof(double));
      memcpy(vec_sphere_exs + offset, rhs.vec_sphere_exs, chunk_size * sizeof(double));
      memcpy(vec_sphere_albs + offset, rhs.vec_sphere_albs, chunk_size * sizeof(double));
      memcpy(vec_sphere_sqscs + offset, rhs.vec_sphere_sqscs, chunk_size * sizeof(double));
      memcpy(vec_sphere_sqabs + offset, rhs.vec_sphere_sqabs, chunk_size * sizeof(double));
      memcpy(vec_sphere_sqexs + offset, rhs.vec_sphere_sqexs, chunk_size * sizeof(double));
      memcpy(vec_fsas + offset, rhs.vec_fsas, chunk_size * sizeof(dcomplex));
      memcpy(vec_qschus + offset, rhs.vec_qschus, chunk_size * sizeof(double));
      memcpy(vec_pschus + offset, rhs.vec_pschus, chunk_size * sizeof(double));
      memcpy(vec_s0mags + offset, rhs.vec_s0mags, chunk_size * sizeof(double));
      memcpy(vec_cosavs + offset, rhs.vec_cosavs, chunk_size * sizeof(double));
      memcpy(vec_raprs + offset, rhs.vec_raprs, chunk_size * sizeof(double));
      memcpy(vec_tqek1 + offset, rhs.vec_tqek1, chunk_size * sizeof(double));
      memcpy(vec_tqsk1 + offset, rhs.vec_tqsk1, chunk_size * sizeof(double));
      memcpy(vec_tqek2 + offset, rhs.vec_tqek2, chunk_size * sizeof(double));
      memcpy(vec_tqsk2 + offset, rhs.vec_tqsk2, chunk_size * sizeof(double));
      // Insert vectors of multiple directions per configuration
      offset = (xi1 - _first_xi) * configurations * ndirs;
      chunk_size = rhs.xi_block_size * configurations * ndirs;
      memcpy(vec_dir_sas11 + offset, rhs.vec_dir_sas11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas21 + offset, rhs.vec_dir_sas21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas12 + offset, rhs.vec_dir_sas12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas22 + offset, rhs.vec_dir_sas22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_muls + 16 * offset, rhs.vec_dir_muls, 16 * chunk_size * sizeof(double));
      memcpy(vec_dir_mulslr + 16 * offset, rhs.vec_dir_mulslr, 16 * chunk_size * sizeof(double));
      // Insert vectors whose sizes depend on wavelengths and directions
      offset = (xi1 - _first_xi) * ndirs;
      chunk_size = rhs.xi_block_size * ndirs;
      memcpy(vec_dir_sat11 + offset, rhs.vec_dir_sat11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sat21 + offset, rhs.vec_dir_sat21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sat12 + offset, rhs.vec_dir_sat12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sat22 + offset, rhs.vec_dir_sat22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_scc1 + offset, rhs.vec_dir_scc1, chunk_size * sizeof(double));
      memcpy(vec_dir_scc2 + offset, rhs.vec_dir_scc2, chunk_size * sizeof(double));
      memcpy(vec_dir_abc1 + offset, rhs.vec_dir_abc1, chunk_size * sizeof(double));
      memcpy(vec_dir_abc2 + offset, rhs.vec_dir_abc2, chunk_size * sizeof(double));
      memcpy(vec_dir_exc1 + offset, rhs.vec_dir_exc1, chunk_size * sizeof(double));
      memcpy(vec_dir_exc2 + offset, rhs.vec_dir_exc2, chunk_size * sizeof(double));
      memcpy(vec_dir_albedc1 + offset, rhs.vec_dir_albedc1, chunk_size * sizeof(double));
      memcpy(vec_dir_albedc2 + offset, rhs.vec_dir_albedc2, chunk_size * sizeof(double));
      memcpy(vec_dir_qscc1 + offset, rhs.vec_dir_qscc1, chunk_size * sizeof(double));
      memcpy(vec_dir_qscc2 + offset, rhs.vec_dir_qscc2, chunk_size * sizeof(double));
      memcpy(vec_dir_qabc1 + offset, rhs.vec_dir_qabc1, chunk_size * sizeof(double));
      memcpy(vec_dir_qabc2 + offset, rhs.vec_dir_qabc2, chunk_size * sizeof(double));
      memcpy(vec_dir_qexc1 + offset, rhs.vec_dir_qexc1, chunk_size * sizeof(double));
      memcpy(vec_dir_qexc2 + offset, rhs.vec_dir_qexc2, chunk_size * sizeof(double));
      memcpy(vec_dir_sccrt1 + offset, rhs.vec_dir_sccrt1, chunk_size * sizeof(double));
      memcpy(vec_dir_sccrt2 + offset, rhs.vec_dir_sccrt2, chunk_size * sizeof(double));
      memcpy(vec_dir_abcrt1 + offset, rhs.vec_dir_abcrt1, chunk_size * sizeof(double));
      memcpy(vec_dir_abcrt2 + offset, rhs.vec_dir_abcrt2, chunk_size * sizeof(double));
      memcpy(vec_dir_excrt1 + offset, rhs.vec_dir_excrt1, chunk_size * sizeof(double));
      memcpy(vec_dir_excrt2 + offset, rhs.vec_dir_excrt2, chunk_size * sizeof(double));
      memcpy(vec_dir_fsac11 + offset, rhs.vec_dir_fsac11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_fsac21 + offset, rhs.vec_dir_fsac21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_fsac12 + offset, rhs.vec_dir_fsac12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_fsac22 + offset, rhs.vec_dir_fsac22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sac11 + offset, rhs.vec_dir_sac11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sac21 + offset, rhs.vec_dir_sac21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sac12 + offset, rhs.vec_dir_sac12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sac22 + offset, rhs.vec_dir_sac22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_qschuc1 + offset, rhs.vec_dir_qschuc1, chunk_size * sizeof(double));
      memcpy(vec_dir_qschuc2 + offset, rhs.vec_dir_qschuc2, chunk_size * sizeof(double));
      memcpy(vec_dir_pschuc1 + offset, rhs.vec_dir_pschuc1, chunk_size * sizeof(double));
      memcpy(vec_dir_pschuc2 + offset, rhs.vec_dir_pschuc2, chunk_size * sizeof(double));
      memcpy(vec_dir_s0magc1 + offset, rhs.vec_dir_s0magc1, chunk_size * sizeof(double));
      memcpy(vec_dir_s0magc2 + offset, rhs.vec_dir_s0magc2, chunk_size * sizeof(double));
      memcpy(vec_dir_cosavc1 + offset, rhs.vec_dir_cosavc1, chunk_size * sizeof(double));
      memcpy(vec_dir_cosavc2 + offset, rhs.vec_dir_cosavc2, chunk_size * sizeof(double));
      memcpy(vec_dir_raprc1 + offset, rhs.vec_dir_raprc1, chunk_size * sizeof(double));
      memcpy(vec_dir_raprc2 + offset, rhs.vec_dir_raprc2, chunk_size * sizeof(double));
      memcpy(vec_dir_flc1 + offset, rhs.vec_dir_flc1, chunk_size * sizeof(double));
      memcpy(vec_dir_flc2 + offset, rhs.vec_dir_flc2, chunk_size * sizeof(double));
      memcpy(vec_dir_frc1 + offset, rhs.vec_dir_frc1, chunk_size * sizeof(double));
      memcpy(vec_dir_frc2 + offset, rhs.vec_dir_frc2, chunk_size * sizeof(double));
      memcpy(vec_dir_fkc1 + offset, rhs.vec_dir_fkc1, chunk_size * sizeof(double));
      memcpy(vec_dir_fkc2 + offset, rhs.vec_dir_fkc2, chunk_size * sizeof(double));
      memcpy(vec_dir_fxc1 + offset, rhs.vec_dir_fxc1, chunk_size * sizeof(double));
      memcpy(vec_dir_fxc2 + offset, rhs.vec_dir_fxc2, chunk_size * sizeof(double));
      memcpy(vec_dir_fyc1 + offset, rhs.vec_dir_fyc1, chunk_size * sizeof(double));
      memcpy(vec_dir_fyc2 + offset, rhs.vec_dir_fyc2, chunk_size * sizeof(double));
      memcpy(vec_dir_fzc1 + offset, rhs.vec_dir_fzc1, chunk_size * sizeof(double));
      memcpy(vec_dir_fzc2 + offset, rhs.vec_dir_fzc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqelc1 + offset, rhs.vec_dir_tqelc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqelc2 + offset, rhs.vec_dir_tqelc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqerc1 + offset, rhs.vec_dir_tqerc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqerc2 + offset, rhs.vec_dir_tqerc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqekc1 + offset, rhs.vec_dir_tqekc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqekc2 + offset, rhs.vec_dir_tqekc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqexc1 + offset, rhs.vec_dir_tqexc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqexc2 + offset, rhs.vec_dir_tqexc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqeyc1 + offset, rhs.vec_dir_tqeyc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqeyc2 + offset, rhs.vec_dir_tqeyc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqezc1 + offset, rhs.vec_dir_tqezc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqezc2 + offset, rhs.vec_dir_tqezc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqslc1 + offset, rhs.vec_dir_tqslc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqslc2 + offset, rhs.vec_dir_tqslc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsrc1 + offset, rhs.vec_dir_tqsrc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsrc2 + offset, rhs.vec_dir_tqsrc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqskc1 + offset, rhs.vec_dir_tqskc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqskc2 + offset, rhs.vec_dir_tqskc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsxc1 + offset, rhs.vec_dir_tqsxc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsxc2 + offset, rhs.vec_dir_tqsxc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsyc1 + offset, rhs.vec_dir_tqsyc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsyc2 + offset, rhs.vec_dir_tqsyc2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqszc1 + offset, rhs.vec_dir_tqszc1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqszc2 + offset, rhs.vec_dir_tqszc2, chunk_size * sizeof(double));
      memcpy(vec_dir_mulc + 16 * offset, rhs.vec_dir_mulc, 16 * chunk_size * sizeof(double));
      memcpy(vec_dir_mulclr + 16 * offset, rhs.vec_dir_mulclr, 16 * chunk_size * sizeof(double));
    }
  }
  return result;
}

int ClusterOutputInfo::write(const std::string &output, const std::string &format) {
  int result = 0;
  if (format.compare("LEGACY") == 0) {
    result = write_legacy(output);
  } else if (format.compare("HDF5") == 0) {
    result = write_hdf5(output);
  } else {
    string message = "Unknown format mode: \"" + format + "\"";
    throw UnrecognizedConfigurationException(message);
  }
  return result;
}

int ClusterOutputInfo::write_hdf5(const std::string &file_name) {
  List<string> *rec_name_list = new List<string>(1);
  List<string> *rec_type_list = new List<string>(1);
  List<void *> *rec_ptr_list = new List<void *>(1);
  string str_type, str_name;
  rec_name_list->set(0, "NSPH");
  rec_type_list->set(0, "INT32_(1)");
  rec_ptr_list->set(0, &nsph);
  rec_name_list->append("LI");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&li);
  rec_name_list->append("LE");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&le);
  rec_name_list->append("LM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&lm);
  rec_name_list->append("MXNDM");
  rec_type_list->append("INT64_(1)");
  rec_ptr_list->append(&mxndm);
  rec_name_list->append("INPOL");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&inpol);
  rec_name_list->append("NPNT");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&npnt);
  rec_name_list->append("NPNTTS");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&npntts);
  rec_name_list->append("IAVM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&iavm);
  rec_name_list->append("ISAM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&isam);
  rec_name_list->append("JWTM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&jwtm);
  rec_name_list->append("VEC_SPH_X");
  rec_type_list->append("FLOAT64_(" + to_string(nsph) + ")");
  rec_ptr_list->append(vec_x_coords);
  rec_name_list->append("VEC_SPH_Y");
  rec_type_list->append("FLOAT64_(" + to_string(nsph) + ")");
  rec_ptr_list->append(vec_y_coords);
  rec_name_list->append("VEC_SPH_Z");
  rec_type_list->append("FLOAT64_(" + to_string(nsph) + ")");
  rec_ptr_list->append(vec_z_coords);
  rec_name_list->append("TH_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&th);
  rec_name_list->append("TH_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thstp);
  rec_name_list->append("TH_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thlst);
  rec_name_list->append("THS_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&ths);
  rec_name_list->append("THS_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thsstp);
  rec_name_list->append("THS_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thslst);
  rec_name_list->append("PH_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&ph);
  rec_name_list->append("PH_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phstp);
  rec_name_list->append("PH_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phlst);
  rec_name_list->append("PHS_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phs);
  rec_name_list->append("PHS_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phsstp);
  rec_name_list->append("PHS_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phslst);
  rec_name_list->append("EXRI");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&exri);
  rec_name_list->append("IDFC");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&idfc);
  rec_name_list->append("XI1");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&_first_xi);
  rec_name_list->append("NXI");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&xi_block_size);
  rec_name_list->append("VEC_JXI");
  rec_type_list->append("INT32_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_jxi);
  rec_name_list->append("VEC_IER");
  rec_type_list->append("INT16_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_ier);
  rec_name_list->append("VEC_VK");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_vk);
  rec_name_list->append("VEC_XI");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_xi);
  rec_name_list->append("NCONF");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&configurations);
  rec_name_list->append("VEC_SPH_SIZES");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_sizes);
  rec_name_list->append("VEC_SPH_REFRI");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_ref_indices);
  rec_name_list->append("VEC_SPH_SCS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_scs);
  rec_name_list->append("VEC_SPH_ABS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_abs);
  rec_name_list->append("VEC_SPH_EXS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_exs);
  rec_name_list->append("VEC_SPH_ALBS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_albs);
  rec_name_list->append("VEC_SPH_SQSCS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_sqscs);
  rec_name_list->append("VEC_SPH_SQABS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_sqabs);
  rec_name_list->append("VEC_SPH_SQEXS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_sphere_sqexs);
  rec_name_list->append("VEC_FSAS");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_fsas);
  rec_name_list->append("VEC_QSCHUS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_qschus);
  rec_name_list->append("VEC_PSCHUS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_pschus);
  rec_name_list->append("VEC_S0MAGS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_s0mags);
  rec_name_list->append("VEC_COSAVS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_cosavs);
  rec_name_list->append("VEC_RAPRS");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_raprs);
  rec_name_list->append("VEC_TQEK1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_tqek1);
  rec_name_list->append("VEC_TQSK1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_tqsk1);
  rec_name_list->append("VEC_TQEK2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_tqek2);
  rec_name_list->append("VEC_TQSK2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * configurations) + ")");
  rec_ptr_list->append(vec_tqsk2);
  rec_name_list->append("VEC_FSAT");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size) + ")");
  rec_ptr_list->append(vec_fsat);
  rec_name_list->append("VEC_QSCHUT");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qschut);
  rec_name_list->append("VEC_PSCHUT");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_pschut);
  rec_name_list->append("VEC_S0MAGT");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_s0magt);
  rec_name_list->append("VEC_SCC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_scc1);
  rec_name_list->append("VEC_SCC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_scc2);
  rec_name_list->append("VEC_ABC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_abc1);
  rec_name_list->append("VEC_ABC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_abc2);
  rec_name_list->append("VEC_EXC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_exc1);
  rec_name_list->append("VEC_EXC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_exc2);
  rec_name_list->append("VEC_ALBEDC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_albedc1);
  rec_name_list->append("VEC_ALBEDC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_albedc2);
  rec_name_list->append("VEC_QSCAMC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qscamc1);
  rec_name_list->append("VEC_QSCAMC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qscamc2);
  rec_name_list->append("VEC_QABSMC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qabsmc1);
  rec_name_list->append("VEC_QABSMC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qabsmc2);
  rec_name_list->append("VEC_QEXTMC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qextmc1);
  rec_name_list->append("VEC_QEXTMC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qextmc2);
  rec_name_list->append("VEC_SCCRT1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_sccrt1);
  rec_name_list->append("VEC_SCCRT2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_sccrt2);
  rec_name_list->append("VEC_ABCRT1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_abcrt1);
  rec_name_list->append("VEC_ABCRT2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_abcrt2);
  rec_name_list->append("VEC_EXCRT1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_excrt1);
  rec_name_list->append("VEC_EXCRT2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_excrt2);
  rec_name_list->append("VEC_FSAC11");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size) + ")");
  rec_ptr_list->append(vec_fsac11);
  rec_name_list->append("VEC_FSAC21");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size) + ")");
  rec_ptr_list->append(vec_fsac21);
  rec_name_list->append("VEC_FSAC22");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size) + ")");
  rec_ptr_list->append(vec_fsac22);
  rec_name_list->append("VEC_FSAC12");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size) + ")");
  rec_ptr_list->append(vec_fsac12);
  rec_name_list->append("VEC_QSCHUC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qschuc1);
  rec_name_list->append("VEC_QSCHUC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_qschuc2);
  rec_name_list->append("VEC_PSCHUC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_pschuc1);
  rec_name_list->append("VEC_PSCHUC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_pschuc2);
  rec_name_list->append("VEC_S0MAGC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_s0magc1);
  rec_name_list->append("VEC_S0MAGC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_s0magc2);
  rec_name_list->append("VEC_COSAVC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_cosavc1);
  rec_name_list->append("VEC_COSAVC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_cosavc2);
  rec_name_list->append("VEC_RAPRC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_raprc1);
  rec_name_list->append("VEC_RAPRC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_raprc2);
  rec_name_list->append("VEC_FKC1");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_fkc1);
  rec_name_list->append("VEC_FKC2");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_fkc2);
  rec_name_list->append("VEC_DIR_SCAN");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs) + ")");
  rec_ptr_list->append(vec_dir_scand);
  rec_name_list->append("VEC_DIR_CFMP");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs) + ")");
  rec_ptr_list->append(vec_dir_cfmp);
  rec_name_list->append("VEC_DIR_SFMP");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs) + ")");
  rec_ptr_list->append(vec_dir_sfmp);
  rec_name_list->append("VEC_DIR_CFSP");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs) + ")");
  rec_ptr_list->append(vec_dir_cfsp);
  rec_name_list->append("VEC_DIR_SFSP");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs) + ")");
  rec_ptr_list->append(vec_dir_sfsp);
  rec_name_list->append("VEC_DIR_UN");
  rec_type_list->append("FLOAT64_(" + to_string(3 * ndirs) + ")");
  rec_ptr_list->append(vec_dir_un);
  rec_name_list->append("VEC_DIR_UNS");
  rec_type_list->append("FLOAT64_(" + to_string(3 * ndirs) + ")");
  rec_ptr_list->append(vec_dir_uns);
  rec_name_list->append("VEC_DIR_SAS11");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sas11);
  rec_name_list->append("VEC_DIR_SAS21");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sas21);
  rec_name_list->append("VEC_DIR_SAS12");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sas12);
  rec_name_list->append("VEC_DIR_SAS22");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sas22);
  rec_name_list->append("VEC_DIR_MULS");
  rec_type_list->append("FLOAT64_(" + to_string(16 * ndirs * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_muls);
  rec_name_list->append("VEC_DIR_MULSLR");
  rec_type_list->append("FLOAT64_(" + to_string(16 * ndirs * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_mulslr);
  rec_name_list->append("VEC_DIR_SAT11");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sat11);
  rec_name_list->append("VEC_DIR_SAT21");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sat21);
  rec_name_list->append("VEC_DIR_SAT12");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sat12);
  rec_name_list->append("VEC_DIR_SAT22");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sat22);
  rec_name_list->append("VEC_DIR_SCC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_scc1);
  rec_name_list->append("VEC_DIR_SCC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_scc2);
  rec_name_list->append("VEC_DIR_ABC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_abc1);
  rec_name_list->append("VEC_DIR_ABC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_abc2);
  rec_name_list->append("VEC_DIR_EXC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_exc1);
  rec_name_list->append("VEC_DIR_EXC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_exc2);
  rec_name_list->append("VEC_DIR_ALBEDC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_albedc1);
  rec_name_list->append("VEC_DIR_ALBEDC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_albedc2);
  rec_name_list->append("VEC_DIR_QSCC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qscc1);
  rec_name_list->append("VEC_DIR_QSCC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qscc2);
  rec_name_list->append("VEC_DIR_QABC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qabc1);
  rec_name_list->append("VEC_DIR_QABC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qabc2);
  rec_name_list->append("VEC_DIR_QEXC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qexc1);
  rec_name_list->append("VEC_DIR_QEXC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qexc2);
  rec_name_list->append("VEC_DIR_SCCRT1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sccrt1);
  rec_name_list->append("VEC_DIR_SCCRT2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sccrt2);
  rec_name_list->append("VEC_DIR_ABCRT1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_abcrt1);
  rec_name_list->append("VEC_DIR_ABCRT2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_abcrt2);
  rec_name_list->append("VEC_DIR_EXCRT1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_excrt1);
  rec_name_list->append("VEC_DIR_EXCRT2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_excrt2);
  rec_name_list->append("VEC_DIR_FSAC11");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fsac11);
  rec_name_list->append("VEC_DIR_FSAC21");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fsac21);
  rec_name_list->append("VEC_DIR_FSAC12");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fsac12);
  rec_name_list->append("VEC_DIR_FSAC22");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fsac22);
  rec_name_list->append("VEC_DIR_SAC11");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sac11);
  rec_name_list->append("VEC_DIR_SAC21");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sac21);
  rec_name_list->append("VEC_DIR_SAC12");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sac12);
  rec_name_list->append("VEC_DIR_SAC22");
  rec_type_list->append("FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_sac22);
  rec_name_list->append("VEC_DIR_QSCHUC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qschuc1);
  rec_name_list->append("VEC_DIR_QSCHUC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_qschuc2);
  rec_name_list->append("VEC_DIR_PSCHUC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_pschuc1);
  rec_name_list->append("VEC_DIR_PSCHUC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_pschuc2);
  rec_name_list->append("VEC_DIR_S0MAGC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_s0magc1);
  rec_name_list->append("VEC_DIR_S0MAGC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_s0magc2);
  rec_name_list->append("VEC_DIR_COSAVC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_cosavc1);
  rec_name_list->append("VEC_DIR_COSAVC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_cosavc2);
  rec_name_list->append("VEC_DIR_RAPRC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_raprc1);
  rec_name_list->append("VEC_DIR_RAPRC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_raprc2);
  rec_name_list->append("VEC_DIR_FLC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_flc1);
  rec_name_list->append("VEC_DIR_FLC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_flc2);
  rec_name_list->append("VEC_DIR_FRC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_frc1);
  rec_name_list->append("VEC_DIR_FRC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_frc2);
  rec_name_list->append("VEC_DIR_FKC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fkc1);
  rec_name_list->append("VEC_DIR_FKC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fkc2);
  rec_name_list->append("VEC_DIR_FXC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fxc1);
  rec_name_list->append("VEC_DIR_FXC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fxc2);
  rec_name_list->append("VEC_DIR_FYC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fyc1);
  rec_name_list->append("VEC_DIR_FYC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fyc2);
  rec_name_list->append("VEC_DIR_FZC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fzc1);
  rec_name_list->append("VEC_DIR_FZC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_fzc2);
  rec_name_list->append("VEC_DIR_TQELC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqelc1);
  rec_name_list->append("VEC_DIR_TQELC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqelc2);
  rec_name_list->append("VEC_DIR_TQERC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqerc1);
  rec_name_list->append("VEC_DIR_TQERC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqerc2);
  rec_name_list->append("VEC_DIR_TQEKC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqekc1);
  rec_name_list->append("VEC_DIR_TQEKC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqekc2);
  rec_name_list->append("VEC_DIR_TQEXC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqexc1);
  rec_name_list->append("VEC_DIR_TQEXC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqexc2);
  rec_name_list->append("VEC_DIR_TQEYC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqeyc1);
  rec_name_list->append("VEC_DIR_TQEYC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqeyc2);
  rec_name_list->append("VEC_DIR_TQEZC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqezc1);
  rec_name_list->append("VEC_DIR_TQEZC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqezc2);
  rec_name_list->append("VEC_DIR_TQSLC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqslc1);
  rec_name_list->append("VEC_DIR_TQSLC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqslc2);
  rec_name_list->append("VEC_DIR_TQSRC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqsrc1);
  rec_name_list->append("VEC_DIR_TQSRC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqsrc2);
  rec_name_list->append("VEC_DIR_TQSKC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqskc1);
  rec_name_list->append("VEC_DIR_TQSKC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqskc2);
  rec_name_list->append("VEC_DIR_TQSXC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqsxc1);
  rec_name_list->append("VEC_DIR_TQSXC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqsxc2);
  rec_name_list->append("VEC_DIR_TQSYC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqsyc1);
  rec_name_list->append("VEC_DIR_TQSYC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqsyc2);
  rec_name_list->append("VEC_DIR_TQSZC1");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqszc1);
  rec_name_list->append("VEC_DIR_TQSZC2");
  rec_type_list->append("FLOAT64_(" + to_string(ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_tqszc2);
  rec_name_list->append("VEC_DIR_MULC");
  rec_type_list->append("FLOAT64_(" + to_string(16 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_mulc);
  rec_name_list->append("VEC_DIR_MULCLR");
  rec_type_list->append("FLOAT64_(" + to_string(16 * ndirs * xi_block_size) + ")");
  rec_ptr_list->append(vec_dir_mulclr);

  // Convert the lists to arrays and write them to HDF5
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
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete schema;
  delete hdf_file;
  return 0;
}

int ClusterOutputInfo::write_legacy(const std::string &output) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  int result = 0;
  FILE *p_outfile = fopen(output.c_str(), "w");
  if (p_outfile != NULL) {
    if (vec_jxi[0] == 1) {
      // Write the preamble of c_OCLU.
      fprintf(p_outfile, " READ(IR,*)NSPH,LI,LE,MXNDM,INPOL,NPNT,NPNTTS,IAVM,ISAM\n");
#ifdef USE_ILP64
      fprintf(
	p_outfile, " %5d%5d%5d%5ld%5d%5d%5d%5d%5d\n",
	nsph, li, le, mxndm, inpol, npnt, npntts, iavm, isam
      );
#else
      fprintf(
	p_outfile, " %5d%5d%5d%5d%5d%5d%5d%5d%5d\n",
	nsph, li, le, mxndm, inpol, npnt, npntts, iavm, isam
      );
#endif // USE_ILP64
      fprintf(p_outfile, " READ(IR,*)RXX(I),RYY(I),RZZ(I)\n");
      for (int ri = 0; ri < nsph; ri++) {
	fprintf(
	  p_outfile, "%17.8lE%17.8lE%17.8lE\n",
	  vec_x_coords[ri], vec_y_coords[ri], vec_z_coords[ri]
	);
      }
      fprintf(p_outfile, " READ(IR,*)TH,THSTP,THLST,THS,THSSTP,THSLST\n");
      fprintf(
	p_outfile, " %10.3lE%10.3lE%10.3lE%10.3lE%10.3lE%10.3lE\n",
	th, thstp, thlst, ths, thsstp, thslst
      );
      fprintf(p_outfile, " READ(IR,*)PH,PHSTP,PHLST,PHS,PHSSTP,PHSLST\n");
      fprintf(
	p_outfile, " %10.3lE%10.3lE%10.3lE%10.3lE%10.3lE%10.3lE\n",
	ph, phstp, phlst, phs, phsstp, phslst
      );
      fprintf(p_outfile, " READ(IR,*)JWTM\n");
      fprintf(p_outfile, " %5d\n", jwtm);
      fprintf(p_outfile, "  READ(ITIN)NSPHT\n");
      fprintf(p_outfile, "  READ(ITIN)(IOG(I),I=1,NSPH)\n");
      fprintf(p_outfile, "  READ(ITIN)EXDC,WP,XIP,IDFC,NXI\n");
      fprintf(p_outfile, "  READ(ITIN)(XIV(I),I=1,NXI)\n");
      fprintf(p_outfile, "  READ(ITIN)NSHL(I),ROS(I)\n");
      fprintf(p_outfile, "  READ(ITIN)(RCF(I,NS),NS=1,NSH)\n \n");
      fprintf(p_outfile, "  REFR. INDEX OF EXTERNAL MEDIUM=%15.7lE\n", exri);
      if (idfc < 0) {
	fprintf(p_outfile, "  VK=%15.7lE, XI IS SCALE FACTOR FOR LENGTHS\n \n", vec_vk[0]);
      }
    }
    // Write the wavelength loop data.
    for (int jxi = 0; jxi < xi_block_size; jxi++) {
      fprintf(p_outfile, "========== JXI =%3d ====================\n", vec_jxi[jxi]);
      if (idfc >= 0) {
	fprintf(p_outfile, "  VK=%15.7lE, XI=%15.7lE\n", vec_vk[jxi], vec_xi[jxi]);
      } else {
	fprintf(p_outfile, "  XI=%15.7lE\n", vec_xi[jxi]);
      }
      if (vec_ier[jxi] == 1) {
	fprintf(p_outfile, "  STOP IN HJV\n");
	result = 1;
	break; // jxi loop
      }
      if (vec_ier[jxi] == 2) {
	fprintf(p_outfile, "  STOP IN DME\n");
	result = 2;
	break; // jxi loop
      }
      double alamb = 2.0 * 3.141592653589793 / vec_vk[jxi];
      if (inpol == 0)
	fprintf(p_outfile, "   LIN\n");
      else
	fprintf(p_outfile, "  CIRC\n");
      if (li != le)
	fprintf(p_outfile, "     SPHERES; LMX=LI\n");
      for (int si = 0; si < configurations; si++) {
	fprintf(p_outfile, "     SPHERE %2d\n", si + 1);
	if (vec_sphere_ref_indices[jxi * configurations + si] == cc0)
	  fprintf(p_outfile, "  SIZE=%15.7lE\n", vec_sphere_sizes[jxi * configurations + si]);
	else
	  fprintf(
		  p_outfile, "  SIZE=%15.7lE, REFRACTIVE INDEX=%15.7lE%15.7lE\n",
		  vec_sphere_sizes[jxi * configurations + si],
		  real(vec_sphere_ref_indices[jxi * configurations + si]),
		  imag(vec_sphere_ref_indices[jxi * configurations + si])
	  );
	fprintf(p_outfile, " ----- SCS ----- ABS ----- EXS ----- ALBEDS --\n");
	fprintf(
		p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
		vec_sphere_scs[jxi * configurations + si],
		vec_sphere_abs[jxi * configurations + si],
		vec_sphere_exs[jxi * configurations + si],
		vec_sphere_albs[jxi * configurations + si]
        );
	fprintf(p_outfile, " ---- SCS/GS -- ABS/GS -- EXS/GS ---\n");
	fprintf(
		p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		vec_sphere_sqscs[jxi * configurations + si],
		vec_sphere_sqabs[jxi * configurations + si],
		vec_sphere_sqexs[jxi * configurations + si]
        );
	fprintf(
		p_outfile, "  FSAS=%15.7lE%15.7lE\n",
		real(vec_fsas[jxi * configurations + si]),
		imag(vec_fsas[jxi * configurations + si])
        );
	// fprintf(
	// 	p_outfile, "INSERTION: CS_SPHERE  %15.7lE%15.7lE%15.7lE%15.7lE\n",
	// 	alamb, vec_sphere_scs[jxi * configurations + si],
	// 	vec_sphere_abs[jxi * configurations + si],
	// 	vec_sphere_exs[jxi * configurations + si]
        // );
	fprintf(
		p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
		vec_qschus[jxi * configurations + si],
		vec_pschus[jxi * configurations + si],
		vec_s0mags[jxi * configurations + si]
	);
	fprintf(
		p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
		vec_cosavs[jxi * configurations + si],
		vec_raprs[jxi * configurations + si]
	);
	fprintf(
		p_outfile, "  IPO= 1, TQEk=%15.7lE, TQSk=%15.7lE\n",
		vec_tqek1[jxi * configurations + si],
		vec_tqsk1[jxi * configurations + si]
	);
	fprintf(
		p_outfile, "  IPO= 2, TQEk=%15.7lE, TQSk=%15.7lE\n",
		vec_tqek2[jxi * configurations + si],
		vec_tqsk2[jxi * configurations + si]
	);
      } // si configuration loop
      fprintf(
	      p_outfile, "  FSAT=%15.7lE%15.7lE\n",
	      real(vec_fsat[jxi]), imag(vec_fsat[jxi])
      );
      fprintf(
	      p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	      vec_qschut[jxi], vec_pschut[jxi], vec_s0magt[jxi]
      );
      fprintf(p_outfile, "     CLUSTER (ENSEMBLE AVERAGE, MODE%2d)\n", iavm);
      // Parallel polarization cluster average section
      if (inpol == 0)
	fprintf(p_outfile, "   LIN -1\n");
      else
	fprintf(p_outfile, "  CIRC -1\n");
      fprintf(p_outfile, " ----- SCC ----- ABC ----- EXC ----- ALBEDC --\n");
      fprintf(
	      p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
	      vec_scc1[jxi], vec_abc1[jxi], vec_exc1[jxi], vec_albedc1[jxi]
      );
      fprintf(p_outfile, " --- SCC/TGS - ABC/TGS - EXC/TGS ---\n");
      fprintf(
	      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	      vec_qscamc1[jxi], vec_qabsmc1[jxi], vec_qextmc1[jxi]
      );
      fprintf(p_outfile, " ---- SCCRT --- ABCRT --- EXCRT ----\n");
      fprintf(
	      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	      vec_sccrt1[jxi], vec_abcrt1[jxi], vec_excrt1[jxi]
      );
      fprintf(
	      p_outfile, "  FSAC(1,1)=%15.7lE%15.7lE   FSAC(2,1)=%15.7lE%15.7lE\n",
	      real(vec_fsac11[jxi]), imag(vec_fsac11[jxi]),
	      real(vec_fsac21[jxi]), imag(vec_fsac21[jxi])
      );
      fprintf(
	      p_outfile, "  RE(FSAC(1,1))/RE(TFSAS)=%15.7lE, IM(FSAC(1,1))/IM(TFSAS)=%15.7lE\n",
	      real(vec_fsac11[jxi]) / real(vec_fsat[jxi]),
	      imag(vec_fsac11[jxi]) / imag(vec_fsat[jxi])
      );
      fprintf(
	      p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	      vec_qschuc1[jxi], vec_pschuc1[jxi], vec_s0magc1[jxi]
      );
      fprintf(
	      p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
	      vec_cosavc1[jxi], vec_raprc1[jxi]
      );
      fprintf(p_outfile, "  Fk=%15.7lE\n", vec_fkc1[jxi]);
      // fprintf(
      // 	      p_outfile, "INSERTION: CSM_CLUSTER  %15.7lE%15.7lE%15.7lE%15.7lE\n",
      // 	      alamb, vec_scc1[jxi], vec_abc1[jxi], vec_exc1[jxi]
      // );
      // Perpendicular polarization cluster average section
      if (inpol == 0)
	fprintf(p_outfile, "   LIN  1\n");
      else
	fprintf(p_outfile, "  CIRC  1\n");
      fprintf(p_outfile, " ----- SCC ----- ABC ----- EXC ----- ALBEDC --\n");
      fprintf(
	      p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
	      vec_scc2[jxi], vec_abc2[jxi], vec_exc2[jxi], vec_albedc2[jxi]
      );
      fprintf(p_outfile, " --- SCC/TGS - ABC/TGS - EXC/TGS ---\n");
      fprintf(
	      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	      vec_qscamc2[jxi], vec_qabsmc2[jxi], vec_qextmc2[jxi]
      );
      fprintf(p_outfile, " ---- SCCRT --- ABCRT --- EXCRT ----\n");
      fprintf(
	      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	      vec_sccrt2[jxi], vec_abcrt2[jxi], vec_excrt2[jxi]
      );
      fprintf(
	      p_outfile, "  FSAC(2,2)=%15.7lE%15.7lE   FSAC(1,2)=%15.7lE%15.7lE\n",
	      real(vec_fsac22[jxi]), imag(vec_fsac22[jxi]),
	      real(vec_fsac12[jxi]), imag(vec_fsac12[jxi])
      );
      fprintf(
	      p_outfile, "  RE(FSAC(2,2))/RE(TFSAS)=%15.7lE, IM(FSAC(2,2))/IM(TFSAS)=%15.7lE\n",
	      real(vec_fsac22[jxi]) / real(vec_fsat[jxi]),
	      imag(vec_fsac22[jxi]) / imag(vec_fsat[jxi])
      );
      fprintf(
	      p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	      vec_qschuc2[jxi], vec_pschuc2[jxi], vec_s0magc2[jxi]
      );
      fprintf(
	      p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
	      vec_cosavc2[jxi], vec_raprc2[jxi]
      );
      fprintf(p_outfile, "  Fk=%15.7lE\n", vec_fkc2[jxi]);
      
      fprintf(
	      p_outfile, "  (RE(FSAC(1,1))-RE(FSAC(2,2)))/RE(FSAC(1,1))=%15.7lE\n",
	      (real(vec_fsac11[jxi]) - real(vec_fsac22[jxi])) / real(vec_fsac11[jxi])
      );
      fprintf(
	      p_outfile, "  (IM(FSAC(1,1))-IM(FSAC(2,2)))/IM(FSAC(1,1))=%15.7lE\n",
	      (imag(vec_fsac11[jxi]) - imag(vec_fsac22[jxi])) / imag(vec_fsac11[jxi])
      );
      // Differential directional loop
      // Loop sorting (outer to inner) is:
      // THETA_INC - PHI_INC - THETA_SCAT - PHI_SCAT
      int dir_index = 0;
      for (int jth = 0; jth < _num_theta; jth++) {
	for (int jph = 0; jph < _num_phi; jph++) {
	  for (int jths = 0; jths < _num_thetas; jths++) {
	    for (int jphs = 0; jphs < _num_phis; jphs++) {
	      fprintf(
		      p_outfile, "********** JTH =%3d, JPH =%3d, JTHS =%3d, JPHS =%3d ********************\n",
		      jth + 1, jph + 1, jths + 1, jphs + 1
	      );
	      fprintf(
		      p_outfile, "  TIDG=%10.3lE, PIDG=%10.3lE, TSDG=%10.3lE, PSDG=%10.3lE\n",
		      th + jth * thstp, ph + jph * phstp,
		      ths + jths * thsstp, phs + jphs * phsstp
	      );
	      fprintf(p_outfile, "  SCAND=%10.3lE\n", vec_dir_scand[dir_index]);
	      fprintf(
		      p_outfile, "  CFMP=%15.7lE, SFMP=%15.7lE\n",
		      vec_dir_cfmp[dir_index], vec_dir_sfmp[dir_index]
	      );
	      fprintf(
		      p_outfile, "  CFSP=%15.7lE, SFSP=%15.7lE\n",
		      vec_dir_cfsp[dir_index], vec_dir_sfsp[dir_index]
	      );
	      if (isam >= 0) {
		fprintf(
			p_outfile, "  UNI=(%12.5lE,%12.5lE,%12.5lE)\n",
			vec_dir_un[3 * dir_index], vec_dir_un[3 * dir_index + 1], vec_dir_un[3 * dir_index + 2]
		);
		fprintf(
			p_outfile, "  UNS=(%12.5lE,%12.5lE,%12.5lE)\n",
			vec_dir_uns[3 * dir_index], vec_dir_uns[3 * dir_index + 1], vec_dir_uns[3 * dir_index + 2]
		);
	      } else
		fprintf(
			p_outfile, "  UN=(%12.5lE,%12.5lE,%12.5lE)\n\n",
			vec_dir_un[3 * dir_index], vec_dir_un[3 * dir_index + 1], vec_dir_un[3 * dir_index + 2]
		);
	      if (inpol == 0)
		fprintf(p_outfile, "   LIN\n");
	      else
		fprintf(p_outfile, "  CIRC\n");
	      if (li != le) fprintf(p_outfile, "     SPHERES; LMX=MIN0(LI,LE)\n");
	      for (int si = 0; si < configurations; si++) {
		int sas_dir_index = jxi * configurations * ndirs + configurations * dir_index + si;
		fprintf(p_outfile, "     SPHERE %2d\n", si + 1);
		fprintf(
			p_outfile, "  SAS(1,1)=%15.7lE%15.7lE, SAS(2,1)=%15.7lE%15.7lE\n",
			real(vec_dir_sas11[sas_dir_index]),
			imag(vec_dir_sas11[sas_dir_index]),
			real(vec_dir_sas21[sas_dir_index]),
			imag(vec_dir_sas21[sas_dir_index])
		);
		fprintf(
			p_outfile, "  SAS(1,2)=%15.7lE%15.7lE, SAS(2,2)=%15.7lE%15.7lE\n",
			real(vec_dir_sas12[sas_dir_index]),
			imag(vec_dir_sas12[sas_dir_index]),
			real(vec_dir_sas22[sas_dir_index]),
			imag(vec_dir_sas22[sas_dir_index])
		);
		fprintf(p_outfile, "  MULS\n");
		for (int i1 = 0; i1 < 4; i1++) {
		  int muls_dir_index = 16 * jxi * ndirs * configurations + 16 * configurations * dir_index + 16 * si + 4 * i1;
		  fprintf(
			  p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
			  vec_dir_muls[muls_dir_index],
			  vec_dir_muls[muls_dir_index + 1],
			  vec_dir_muls[muls_dir_index + 2],
			  vec_dir_muls[muls_dir_index + 3]
		  );
		} // i1 loop
		fprintf(p_outfile, "  MULSLR\n");
		for (int i1 = 0; i1 < 4; i1++) {
		  int muls_dir_index = 16 * jxi * ndirs * configurations + 16 * configurations * dir_index + 16 * si + 4 * i1;
		  fprintf(
			  p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
			  vec_dir_mulslr[muls_dir_index],
			  vec_dir_mulslr[muls_dir_index + 1],
			  vec_dir_mulslr[muls_dir_index + 2],
			  vec_dir_mulslr[muls_dir_index + 3]
		  );
		} // i1 loop
	      } // si loop
	      int sat_dir_index = jxi * ndirs + dir_index;
	      fprintf(
		      p_outfile, "  SAT(1,1)=%15.7lE%15.7lE, SAT(2,1)=%15.7lE%15.7lE\n",
		      real(vec_dir_sat11[sat_dir_index]),
		      imag(vec_dir_sat11[sat_dir_index]),
		      real(vec_dir_sat21[sat_dir_index]),
		      imag(vec_dir_sat21[sat_dir_index])
	      );
	      fprintf(
		      p_outfile, "  SAT(1,2)=%15.7lE%15.7lE, SAT(2,2)=%15.7lE%15.7lE\n",
		      real(vec_dir_sat12[sat_dir_index]),
		      imag(vec_dir_sat12[sat_dir_index]),
		      real(vec_dir_sat22[sat_dir_index]),
		      imag(vec_dir_sat22[sat_dir_index])
	      );
	      bool goto190 = isam >= 0 && (jths > 0 || jphs > 0);
	      fprintf(p_outfile, "     CLUSTER\n");
	      // Parallel polarization cluster section
	      if (inpol == 0)
		fprintf(p_outfile, "   LIN -1\n");
	      else
		fprintf(p_outfile, "  CIRC -1\n");

	      fprintf(p_outfile, " ----- SCC ----- ABC ----- EXC ----- ALBEDC --\n");
	      fprintf(
		      p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
		      vec_dir_scc1[sat_dir_index], vec_dir_abc1[sat_dir_index],
		      vec_dir_exc1[sat_dir_index], vec_dir_albedc1[sat_dir_index]
		      );
	      fprintf(p_outfile, " --- SCC/TGS - ABC/TGS - EXC/TGS ---\n");
	      fprintf(
		      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		      vec_dir_qscc1[sat_dir_index],
		      vec_dir_qabc1[sat_dir_index],
		      vec_dir_qexc1[sat_dir_index]
	      );
	      fprintf(p_outfile, " ---- SCCRT --- ABCRT --- EXCRT ----\n");
	      fprintf(
		      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		      vec_dir_sccrt1[sat_dir_index],
		      vec_dir_abcrt1[sat_dir_index],
		      vec_dir_excrt1[sat_dir_index]
	      );
	      fprintf(
		      p_outfile, "  FSAC(1,1)=%15.7lE%15.7lE   FSAC(2,1)=%15.7lE%15.7lE\n",
		      real(vec_dir_fsac11[sat_dir_index]),
		      imag(vec_dir_fsac11[sat_dir_index]),
		      real(vec_dir_fsac21[sat_dir_index]),
		      imag(vec_dir_fsac21[sat_dir_index])
	      );
	      fprintf(
		      p_outfile, "   SAC(1,1)=%15.7lE%15.7lE    SAC(2,1)=%15.7lE%15.7lE\n",
		      real(vec_dir_sac11[sat_dir_index]),
		      imag(vec_dir_sac11[sat_dir_index]),
		      real(vec_dir_sac21[sat_dir_index]),
		      imag(vec_dir_sac21[sat_dir_index])
	      );
	      fprintf(
		      p_outfile, "  RE(FSAC(1,1))/RE(TFSAS)=%15.7lE, IM(FSAC(1,1))/IM(TFSAS)=%15.7lE\n",
		      real(vec_dir_fsac11[sat_dir_index]) / real(vec_fsat[jxi]),
		      imag(vec_dir_fsac11[sat_dir_index]) / imag(vec_fsat[jxi])
	      );
	      fprintf(
		      p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
		      vec_dir_qschuc1[sat_dir_index],
		      vec_dir_pschuc1[sat_dir_index],
		      vec_dir_s0magc1[sat_dir_index]
	      );
	      // fprintf(
	      // 	      p_outfile, "INSERTION: CS1_CLUSTER  %13.5le%10.3le%10.3le%15.7le%15.7le%15.7le\n",
	      // 	      alamb, th + jth * thstp, ths + jths * thsstp,
	      // 	      vec_dir_scc1[sat_dir_index],
	      // 	      vec_dir_abc1[sat_dir_index],
	      // 	      vec_dir_exc1[sat_dir_index]
	      // );
	      if (!goto190) {
		fprintf(
			p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
			vec_dir_cosavc1[sat_dir_index],
			vec_dir_raprc1[sat_dir_index]
		);
		fprintf(
			p_outfile, "  Fl=%15.7lE, Fr=%15.7lE, Fk=%15.7lE\n",
			vec_dir_flc1[sat_dir_index],
			vec_dir_frc1[sat_dir_index],
			vec_dir_fkc1[sat_dir_index]
		);
		fprintf(
			p_outfile, "  Fx=%15.7lE, Fy=%15.7lE, Fz=%15.7lE\n",
			vec_dir_fxc1[sat_dir_index],
			vec_dir_fyc1[sat_dir_index],
			vec_dir_fzc1[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQEl=%15.7lE,  TQEr=%15.7lE,  TQEk=%15.7lE\n",
			vec_dir_tqelc1[sat_dir_index],
			vec_dir_tqerc1[sat_dir_index],
			vec_dir_tqekc1[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQSl=%15.7lE,  TQSr=%15.7lE,  TQSk=%15.7lE\n",
			vec_dir_tqslc1[sat_dir_index],
			vec_dir_tqsrc1[sat_dir_index],
			vec_dir_tqskc1[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQEx=%15.7lE,  TQEy=%15.7lE,  TQEz=%15.7lE\n",
			vec_dir_tqexc1[sat_dir_index],
			vec_dir_tqeyc1[sat_dir_index],
			vec_dir_tqezc1[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQSx=%15.7lE,  TQSy=%15.7lE,  TQSz=%15.7lE\n",
			vec_dir_tqsxc1[sat_dir_index],
			vec_dir_tqsyc1[sat_dir_index],
			vec_dir_tqszc1[sat_dir_index]
		);
	      } // end goto190 switch
	      // Perpendicular polarization cluster section
	      if (inpol == 0)
		fprintf(p_outfile, "   LIN  1\n");
	      else
		fprintf(p_outfile, "  CIRC  1\n");

	      fprintf(p_outfile, " ----- SCC ----- ABC ----- EXC ----- ALBEDC --\n");
	      fprintf(
		      p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
		      vec_dir_scc2[sat_dir_index], vec_dir_abc2[sat_dir_index],
		      vec_dir_exc2[sat_dir_index], vec_dir_albedc2[sat_dir_index]
		      );
	      fprintf(p_outfile, " --- SCC/TGS - ABC/TGS - EXC/TGS ---\n");
	      fprintf(
		      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		      vec_dir_qscc2[sat_dir_index],
		      vec_dir_qabc2[sat_dir_index],
		      vec_dir_qexc2[sat_dir_index]
	      );
	      fprintf(p_outfile, " ---- SCCRT --- ABCRT --- EXCRT ----\n");
	      fprintf(
		      p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		      vec_dir_sccrt2[sat_dir_index],
		      vec_dir_abcrt2[sat_dir_index],
		      vec_dir_excrt2[sat_dir_index]
	      );
	      fprintf(
		      p_outfile, "  FSAC(2,2)=%15.7lE%15.7lE   FSAC(1,2)=%15.7lE%15.7lE\n",
		      real(vec_dir_fsac22[sat_dir_index]),
		      imag(vec_dir_fsac22[sat_dir_index]),
		      real(vec_dir_fsac12[sat_dir_index]),
		      imag(vec_dir_fsac12[sat_dir_index])
	      );
	      fprintf(
		      p_outfile, "   SAC(2,2)=%15.7lE%15.7lE    SAC(1,2)=%15.7lE%15.7lE\n",
		      real(vec_dir_sac22[sat_dir_index]),
		      imag(vec_dir_sac22[sat_dir_index]),
		      real(vec_dir_sac12[sat_dir_index]),
		      imag(vec_dir_sac12[sat_dir_index])
	      );
	      fprintf(
		      p_outfile, "  RE(FSAC(2,2))/RE(TFSAS)=%15.7lE, IM(FSAC(2,2))/IM(TFSAS)=%15.7lE\n",
		      real(vec_dir_fsac22[sat_dir_index]) / real(vec_fsat[jxi]),
		      imag(vec_dir_fsac22[sat_dir_index]) / imag(vec_fsat[jxi])
	      );
	      fprintf(
		      p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
		      vec_dir_qschuc2[sat_dir_index],
		      vec_dir_pschuc2[sat_dir_index],
		      vec_dir_s0magc2[sat_dir_index]
	      );
	      // fprintf(
	      // 	      p_outfile, "INSERTION: CS2_CLUSTER  %13.5le%10.3le%10.3le%15.7le%15.7le%15.7le\n",
	      // 	      alamb, th + jth * thstp, ths + jths * thsstp,
	      // 	      vec_dir_scc2[sat_dir_index],
	      // 	      vec_dir_abc2[sat_dir_index],
	      // 	      vec_dir_exc2[sat_dir_index]
	      // );
	      if (!goto190) {
		fprintf(
			p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
			vec_dir_cosavc2[sat_dir_index],
			vec_dir_raprc2[sat_dir_index]
		);
		fprintf(
			p_outfile, "  Fl=%15.7lE, Fr=%15.7lE, Fk=%15.7lE\n",
			vec_dir_flc2[sat_dir_index],
			vec_dir_frc2[sat_dir_index],
			vec_dir_fkc2[sat_dir_index]
		);
		fprintf(
			p_outfile, "  Fx=%15.7lE, Fy=%15.7lE, Fz=%15.7lE\n",
			vec_dir_fxc2[sat_dir_index],
			vec_dir_fyc2[sat_dir_index],
			vec_dir_fzc2[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQEl=%15.7lE,  TQEr=%15.7lE,  TQEk=%15.7lE\n",
			vec_dir_tqelc2[sat_dir_index],
			vec_dir_tqerc2[sat_dir_index],
			vec_dir_tqekc2[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQSl=%15.7lE,  TQSr=%15.7lE,  TQSk=%15.7lE\n",
			vec_dir_tqslc2[sat_dir_index],
			vec_dir_tqsrc2[sat_dir_index],
			vec_dir_tqskc2[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQEx=%15.7lE,  TQEy=%15.7lE,  TQEz=%15.7lE\n",
			vec_dir_tqexc2[sat_dir_index],
			vec_dir_tqeyc2[sat_dir_index],
			vec_dir_tqezc2[sat_dir_index]
		);
		fprintf(
			p_outfile, "   TQSx=%15.7lE,  TQSy=%15.7lE,  TQSz=%15.7lE\n",
			vec_dir_tqsxc2[sat_dir_index],
			vec_dir_tqsyc2[sat_dir_index],
			vec_dir_tqszc2[sat_dir_index]
		);
	      } // end goto190 switch
	      fprintf(
		      p_outfile, "  (RE(FSAC(1,1))-RE(FSAC(2,2)))/RE(FSAC(1,1))=%15.7lE\n",
		      (real(vec_dir_fsac11[sat_dir_index]) - real(vec_dir_fsac22[sat_dir_index])) / real(vec_dir_fsac11[sat_dir_index])
	      );
	      fprintf(
		      p_outfile, "  (IM(FSAC(1,1))-IM(FSAC(2,2)))/IM(FSAC(1,1))=%15.7lE\n",
		      (imag(vec_dir_fsac11[sat_dir_index]) - imag(vec_dir_fsac22[sat_dir_index])) / imag(vec_dir_fsac11[sat_dir_index])
	      );
	      fprintf(p_outfile, "  MULC\n");
	      for (int i = 0; i < 4; i++) {
		int mulc_dir_index = 16 * jxi * ndirs + 16 * dir_index + 4 * i;
		fprintf(
			p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
			vec_dir_mulc[mulc_dir_index],
			vec_dir_mulc[mulc_dir_index + 1],
			vec_dir_mulc[mulc_dir_index + 2],
			vec_dir_mulc[mulc_dir_index + 3]
		);
	      } // i mulc loop
	      fprintf(p_outfile, "  MULCLR\n");
	      for (int i = 0; i < 4; i++) {
		int mulc_dir_index = 16 * jxi * ndirs + 16 * dir_index + 4 * i;
		fprintf(
			p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
			vec_dir_mulclr[mulc_dir_index],
			vec_dir_mulclr[mulc_dir_index + 1],
			vec_dir_mulclr[mulc_dir_index + 2],
			vec_dir_mulclr[mulc_dir_index + 3]
		);
	      } // i mulclr loop
	      if (iavm != 0) {
		fprintf(p_outfile, "     CLUSTER (ENSEMBLE AVERAGE, MODE%2d)\n", iavm);
		if (inpol == 0)
		  fprintf(p_outfile, "   LIN\n");
		else
		  fprintf(p_outfile, "  CIRC\n");
		fprintf(p_outfile, "  MULC\n");
		for (int i = 0; i < 4; i++) {
		  int mulc_dir_index = 16 * jxi + 4 * i;
		  fprintf(
			  p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
			  vec_dir_mulc[mulc_dir_index],
			  vec_dir_mulc[mulc_dir_index + 1],
			  vec_dir_mulc[mulc_dir_index + 2],
			  vec_dir_mulc[mulc_dir_index + 3]
		  );
		} // i mulc loop
		fprintf(p_outfile, "  MULCLR\n");
		for (int i = 0; i < 4; i++) {
		  int mulc_dir_index = 16 * jxi + 4 * i;
		  fprintf(
			  p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
			  vec_dir_mulclr[mulc_dir_index],
			  vec_dir_mulclr[mulc_dir_index + 1],
			  vec_dir_mulclr[mulc_dir_index + 2],
			  vec_dir_mulclr[mulc_dir_index + 3]
		  );
		}
	      } // end of if (iavm != 0) switch
	      dir_index++;
	    } // jphs loop
	  } // jths loop
	} // jph loop
      } // jth loop
    } // jxi wavelength loop
    fclose(p_outfile);
  } else {
    result = -1;
  }
  return result;
}

#ifdef MPI_VERSION
int ClusterOutputInfo::mpireceive(const mixMPI *mpidata, int pid) {
  int result = 0;
  int flag;
  int chk_nsph, chk_inpol, chk_iavm, chk_isam, chk_num_theta, chk_num_thetas;
  int chk_num_phi, chk_num_phis, chk_ndirs, chk_idfc, chk_configs;
  double chk_exri;
  MPI_Recv(&flag, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Proceed with the rest _only if_ skip_flag==0, else nothing is to be received
  if (flag == 0) {
    MPI_Recv(&chk_nsph, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_inpol, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_iavm, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_isam, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_theta, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_thetas, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_phi, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_phis, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_ndirs, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_exri, 1, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_idfc, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_configs, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    result += (chk_nsph == nsph) ? 0 : 1;
    result += (chk_inpol == inpol) ? 0 : 1;
    result += (chk_iavm == iavm) ? 0 : 1;
    result += (chk_isam == isam) ? 0 : 1;
    result += (chk_num_theta == _num_theta) ? 0 : 1;
    result += (chk_num_thetas == _num_thetas) ? 0 : 1;
    result += (chk_num_phi == _num_phi) ? 0 : 1;
    result += (chk_num_phis == _num_phis) ? 0 : 1;
    result += (chk_ndirs == ndirs) ? 0 : 1;
    result += (chk_exri == exri) ? 0 : 1;
    result += (chk_idfc == idfc) ? 0 : 1;
    result += (chk_configs == configurations) ? 0 : 1;
    if (result == 0) {
      int xi1, offset, chunk_size;
      MPI_Send(&result, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD);
      MPI_Recv(&xi1, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Receive vectors of single values per scale
      offset = xi1 - _first_xi;
      MPI_Recv(vec_jxi + offset, chunk_size, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_ier + offset, chunk_size, MPI_SHORT, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_vk + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_xi + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsat + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschut + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschut + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0magt + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_albedc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_albedc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qscamc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qscamc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qabsmc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qabsmc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qextmc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qextmc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sccrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sccrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abcrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abcrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_excrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_excrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsac11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsac21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsac22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsac12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschuc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschuc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschuc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschuc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0magc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0magc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_cosavc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_cosavc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_raprc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_raprc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fkc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fkc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors of multiple configuration values per scale
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * configurations;
      MPI_Recv(vec_sphere_sizes + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_ref_indices + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_scs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_abs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_exs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_albs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_sqscs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_sqabs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_sqexs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsas + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschus + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschus + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0mags + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_cosavs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_raprs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqek1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqsk1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqek2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqsk2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors whose sizes depend on directions and configurations.
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * ndirs * configurations;
      MPI_Recv(vec_dir_sas11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_muls + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_mulslr + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors whose sizes depend on directions and scales.
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * ndirs;
      MPI_Recv(vec_dir_sat11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sat21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sat12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sat22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_scc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_scc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_abc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_abc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_exc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_exc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_albedc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_albedc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qscc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qscc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qabc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qabc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qexc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qexc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sccrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sccrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_abcrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_abcrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_excrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_excrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsac11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsac21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsac12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsac22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sac11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sac21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sac12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sac22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qschuc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qschuc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_pschuc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_pschuc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_s0magc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_s0magc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_cosavc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_cosavc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_raprc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_raprc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_flc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_flc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_frc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_frc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fkc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fkc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fxc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fxc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fyc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fyc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fzc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fzc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqelc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqelc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqerc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqerc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqekc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqekc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqexc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqexc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqeyc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqeyc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqezc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqezc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqslc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqslc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsrc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsrc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqskc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqskc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsxc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsxc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsyc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsyc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqszc1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqszc2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_mulc + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_mulclr + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Send(&result, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD);
    }
  }
  return result;
}

int ClusterOutputInfo::mpisend(const mixMPI *mpidata) {
  int result = 0;
  int chunk_size;
  if (_skip_flag == 1) {
    // tell the receiver we are not sending anything
    int flag = 1;
    MPI_Send(&flag, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
  }
  else {
    // tell the receiver we are sending actual stuff
    int flag = 0;
    MPI_Send(&flag, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    // Send output metadata for configuration cross-check
    MPI_Send(&nsph, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&inpol, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&iavm, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&isam, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_theta, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_thetas, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_phi, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_phis, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&ndirs, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&exri, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&idfc, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&configurations, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    // Wait for process 0 to cross-check the configuration
    MPI_Recv(&result, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (result == 0) {
      // Process 0 confirmed the consistency of configuration. Send the data.
      // Send vectors of single values per scale
      MPI_Send(&_first_xi, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(&xi_block_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_jxi, xi_block_size, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_ier, xi_block_size, MPI_SHORT, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_vk, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_xi, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsat, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschut, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschut, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0magt, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_albedc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_albedc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qscamc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qscamc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qabsmc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qabsmc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qextmc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qextmc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sccrt1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sccrt2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abcrt1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abcrt2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_excrt1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_excrt2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsac11, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsac21, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsac22, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsac12, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschuc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschuc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschuc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschuc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0magc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0magc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_cosavc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_cosavc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_raprc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_raprc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fkc1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fkc2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);

      // Send vectors of multiple configuration values per scale
      chunk_size = xi_block_size * configurations;
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_sizes, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_ref_indices, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_scs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_abs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_exs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_albs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_sqscs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_sqabs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_sqexs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsas, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschus, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschus, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0mags, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_cosavs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_raprs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqek1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqsk1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqek2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqsk2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);

      // Send vectors whose sizes depend on directions and configurations.
      chunk_size = ndirs * configurations * xi_block_size;
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_muls, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_mulslr, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);

      // Send vectors whose sizes depend on directions and scales.
      chunk_size = xi_block_size * ndirs;
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sat11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sat21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sat12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sat22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_scc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_scc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_abc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_abc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_exc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_exc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_albedc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_albedc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qscc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qscc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qabc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qabc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qexc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qexc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sccrt1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sccrt2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_abcrt1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_abcrt2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_excrt1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_excrt2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsac11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsac21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsac12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsac22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sac11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sac21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sac12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sac22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qschuc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qschuc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_pschuc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_pschuc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_s0magc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_s0magc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_cosavc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_cosavc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_raprc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_raprc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_flc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_flc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_frc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_frc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fkc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fkc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fxc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fxc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fyc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fyc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fzc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fzc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqelc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqelc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqerc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqerc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqekc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqekc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqexc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqexc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqeyc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqeyc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqezc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqezc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqslc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqslc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsrc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsrc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqskc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqskc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsxc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsxc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsyc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsyc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqszc1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqszc2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_mulc, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_mulclr, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    }
  }
  return result;
}
#endif //MPI_VERSION
// >>> END OF ClusterOutputInfo CLASS IMPLEMENTATION <<<

// >>> InclusionOutputInfo CLASS IMPLEMENTATION <<<
InclusionOutputInfo::InclusionOutputInfo(
  ScattererConfiguration *sc, GeometryConfiguration *gc,
  const mixMPI *mpidata, int first_xi, int xi_length
) {
  _skip_flag = 0;
  nsph = gc->number_of_spheres;
  li = gc->li;
  le = gc->le;
  lm = gc->l_max;
  mxndm = gc->mxndm;
  inpol = gc->in_pol;
  npnt = gc->npnt;
  npntts = gc->npntts;
  iavm = gc->iavm;
  isam = gc->isam;
  jwtm = gc->jwtm;
  // Get spherical constituent coordinates
  vec_x_coords = new double[nsph];
  vec_y_coords = new double[nsph];
  vec_z_coords = new double[nsph];
  for (int nsi = 0; nsi < nsph; nsi++) {
    vec_x_coords[nsi] = gc->get_sph_x(nsi);
    vec_y_coords[nsi] = gc->get_sph_y(nsi);
    vec_z_coords[nsi] = gc->get_sph_z(nsi);
  }
  // Get directional information
  th = gc->in_theta_start;
  thstp = gc->in_theta_step;
  thlst = gc->in_theta_end;
  ths = gc->sc_theta_start;
  thsstp = gc->sc_theta_step;
  thslst = gc->sc_theta_end;
  _num_theta = (thstp == 0.0) ? 1 : 1 + int((thlst - th) / thstp);
  _num_thetas = (thsstp == 0.0) ? 1 : 1 + int((thslst - ths) / thsstp);
  ph = gc->in_phi_start;
  phstp = gc->in_phi_step;
  phlst = gc->in_phi_end;
  phs = gc->sc_phi_start;
  phsstp = gc->sc_phi_step;
  phslst = gc->sc_phi_end;
  _num_phi = (phstp == 0.0) ? 1 : 1 + int((phlst - ph) / phstp);
  _num_phis = (phsstp == 0.0) ? 1 : 1 + int((phslst - phs) / phsstp);
  ndirs = _num_theta * _num_thetas * _num_phi * _num_phis;
  // Get scattering problem configuration
  double exdc = sc->exdc;
  exri = sqrt(exdc);
  idfc = sc->idfc;
  nxi = sc->number_of_scales;
  _first_xi = first_xi;
  xi_block_size = (xi_length == 0) ? nxi - (_first_xi - 1) : xi_length;
  vec_jxi = new int[xi_block_size]();
  vec_ier = new short[xi_block_size]();
  vec_vk = new double[xi_block_size]();
  vec_xi = new double[xi_block_size]();
  configurations = sc->configurations;
  vec_sphere_sizes = new double[xi_block_size * (configurations + 1)]();
  vec_sphere_ref_indices = new dcomplex[xi_block_size * (configurations + 1)]();
  vec_scs1 = new double[xi_block_size]();
  vec_scs2 = new double[xi_block_size]();
  vec_abs1 = new double[xi_block_size]();
  vec_abs2 = new double[xi_block_size]();
  vec_exs1 = new double[xi_block_size]();
  vec_exs2 = new double[xi_block_size]();
  vec_albeds1 = new double[xi_block_size]();
  vec_albeds2 = new double[xi_block_size]();
  vec_scsrt1 = new double[xi_block_size]();
  vec_scsrt2 = new double[xi_block_size]();
  vec_absrt1 = new double[xi_block_size]();
  vec_absrt2 = new double[xi_block_size]();
  vec_exsrt1 = new double[xi_block_size]();
  vec_exsrt2 = new double[xi_block_size]();
  vec_qschu1 = new double[xi_block_size]();
  vec_qschu2 = new double[xi_block_size]();
  vec_pschu1 = new double[xi_block_size]();
  vec_pschu2 = new double[xi_block_size]();
  vec_s0mag1 = new double[xi_block_size]();
  vec_s0mag2 = new double[xi_block_size]();
  vec_cosav1 = new double[xi_block_size]();
  vec_cosav2 = new double[xi_block_size]();
  vec_raprs1 = new double[xi_block_size]();
  vec_raprs2 = new double[xi_block_size]();
  vec_fk1 = new double[xi_block_size]();
  vec_fk2 = new double[xi_block_size]();
  vec_fsas11 = new dcomplex[xi_block_size]();
  vec_fsas21 = new dcomplex[xi_block_size]();
  vec_fsas22 = new dcomplex[xi_block_size]();
  vec_fsas12 = new dcomplex[xi_block_size]();
  vec_dir_tidg = new double[_num_theta];
  vec_dir_pidg = new double[_num_phi];
  vec_dir_tsdg = new double[_num_thetas];
  vec_dir_psdg = new double[_num_phis];
  // Initialize directions (they are scale-independent)
  double cti = th, cpi = ph, cts = ths, cps = phs;
  for (int di = 0; di < _num_theta; di++) {
    vec_dir_tidg[di] = cti;
    cti += thstp;
  }
  for (int di = 0; di < _num_thetas; di++) {
    vec_dir_tsdg[di] = cts;
    cts += thsstp;
  }
  for (int di = 0; di < _num_phi; di++) {
    vec_dir_pidg[di] = cpi;
    cpi += phstp;
  }
  for (int di = 0; di < _num_phis; di++) {
    vec_dir_psdg[di] = cps;
    cps += phsstp;
  }
  vec_dir_scand = new double[ndirs]();
  vec_dir_cfmp = new double[ndirs]();
  vec_dir_sfmp = new double[ndirs]();
  vec_dir_cfsp = new double[ndirs]();
  vec_dir_sfsp = new double[ndirs]();
  vec_dir_un = new double[3 * ndirs]();
  vec_dir_uns = new double[3 * ndirs]();
  vec_dir_scs1 = new double[ndirs * xi_block_size]();
  vec_dir_scs2 = new double[ndirs * xi_block_size]();
  vec_dir_abs1 = new double[ndirs * xi_block_size]();
  vec_dir_abs2 = new double[ndirs * xi_block_size]();
  vec_dir_exs1 = new double[ndirs * xi_block_size]();
  vec_dir_exs2 = new double[ndirs * xi_block_size]();
  vec_dir_albeds1 = new double[ndirs * xi_block_size]();
  vec_dir_albeds2 = new double[ndirs * xi_block_size]();
  vec_dir_scsrt1 = new double[ndirs * xi_block_size]();
  vec_dir_scsrt2 = new double[ndirs * xi_block_size]();
  vec_dir_absrt1 = new double[ndirs * xi_block_size]();
  vec_dir_absrt2 = new double[ndirs * xi_block_size]();
  vec_dir_exsrt1 = new double[ndirs * xi_block_size]();
  vec_dir_exsrt2 = new double[ndirs * xi_block_size]();
  vec_dir_fsas11 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_fsas21 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_fsas12 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_fsas22 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sas11 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sas21 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sas12 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_sas22 = new dcomplex[ndirs * xi_block_size]();
  vec_dir_qschu1 = new double[ndirs * xi_block_size]();
  vec_dir_qschu2 = new double[ndirs * xi_block_size]();
  vec_dir_pschu1 = new double[ndirs * xi_block_size]();
  vec_dir_pschu2 = new double[ndirs * xi_block_size]();
  vec_dir_s0mag1 = new double[ndirs * xi_block_size]();
  vec_dir_s0mag2 = new double[ndirs * xi_block_size]();
  vec_dir_cosav1 = new double[ndirs * xi_block_size]();
  vec_dir_cosav2 = new double[ndirs * xi_block_size]();
  vec_dir_rapr1 = new double[ndirs * xi_block_size]();
  vec_dir_rapr2 = new double[ndirs * xi_block_size]();
  vec_dir_fl1 = new double[ndirs * xi_block_size]();
  vec_dir_fl2 = new double[ndirs * xi_block_size]();
  vec_dir_fr1 = new double[ndirs * xi_block_size]();
  vec_dir_fr2 = new double[ndirs * xi_block_size]();
  vec_dir_fk1 = new double[ndirs * xi_block_size]();
  vec_dir_fk2 = new double[ndirs * xi_block_size]();
  vec_dir_fx1 = new double[ndirs * xi_block_size]();
  vec_dir_fx2 = new double[ndirs * xi_block_size]();
  vec_dir_fy1 = new double[ndirs * xi_block_size]();
  vec_dir_fy2 = new double[ndirs * xi_block_size]();
  vec_dir_fz1 = new double[ndirs * xi_block_size]();
  vec_dir_fz2 = new double[ndirs * xi_block_size]();
  vec_dir_tqel1 = new double[ndirs * xi_block_size]();
  vec_dir_tqel2 = new double[ndirs * xi_block_size]();
  vec_dir_tqer1 = new double[ndirs * xi_block_size]();
  vec_dir_tqer2 = new double[ndirs * xi_block_size]();
  vec_dir_tqek1 = new double[ndirs * xi_block_size]();
  vec_dir_tqek2 = new double[ndirs * xi_block_size]();
  vec_dir_tqex1 = new double[ndirs * xi_block_size]();
  vec_dir_tqex2 = new double[ndirs * xi_block_size]();
  vec_dir_tqey1 = new double[ndirs * xi_block_size]();
  vec_dir_tqey2 = new double[ndirs * xi_block_size]();
  vec_dir_tqez1 = new double[ndirs * xi_block_size]();
  vec_dir_tqez2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsl1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsl2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsr1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsr2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsk1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsk2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsx1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsx2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsy1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsy2 = new double[ndirs * xi_block_size]();
  vec_dir_tqsz1 = new double[ndirs * xi_block_size]();
  vec_dir_tqsz2 = new double[ndirs * xi_block_size]();
  vec_dir_mull = new double[16 * ndirs * xi_block_size]();
  vec_dir_mulllr = new double[16 * ndirs * xi_block_size]();
}

InclusionOutputInfo::InclusionOutputInfo(const std::string &hdf5_name) {
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(hdf5_name, flags);
  herr_t status = hdf_file->get_status();
  string str_name, str_type;
  _skip_flag = 0;
  if (status == 0) {
    status = hdf_file->read("NSPH", "INT32_(1)", &nsph);
    status = hdf_file->read("LI", "INT32_(1)", &li);
    status = hdf_file->read("LE", "INT32_(1)", &le);
    status = hdf_file->read("LM", "INT32_(1)", &lm);
    long tmp;
    status = hdf_file->read("MXNDM", "INT64_(1)", &tmp);
    mxndm = (np_int)tmp;
    status = hdf_file->read("INPOL", "INT32_(1)", &inpol);
    status = hdf_file->read("NPNT", "INT32_(1)", &npnt);
    status = hdf_file->read("NPNTTS", "INT32_(1)", &npntts);
    status = hdf_file->read("IAVM", "INT32_(1)", &iavm);
    status = hdf_file->read("ISAM", "INT32_(1)", &isam);
    status = hdf_file->read("JWTM", "INT32_(1)", &jwtm);
    str_type = "FLOAT64_(" + to_string(nsph) + ")";
    vec_x_coords = new double[nsph];
    vec_y_coords = new double[nsph];
    vec_z_coords = new double[nsph];
    status = hdf_file->read("VEC_SPH_X", str_type, vec_x_coords);
    status = hdf_file->read("VEC_SPH_Y", str_type, vec_y_coords);
    status = hdf_file->read("VEC_SPH_Z", str_type, vec_z_coords);
    status = hdf_file->read("TH_START", "FLOAT64_(1)", &th);
    status = hdf_file->read("TH_STEP", "FLOAT64_(1)", &thstp);
    status = hdf_file->read("TH_END", "FLOAT64_(1)", &thlst);
    _num_theta = (thstp == 0.0) ? 1 : 1 + int((thlst - th) / thstp);
    status = hdf_file->read("THS_START", "FLOAT64_(1)", &ths);
    status = hdf_file->read("THS_STEP", "FLOAT64_(1)", &thsstp);
    status = hdf_file->read("THS_END", "FLOAT64_(1)", &thslst);
    _num_thetas = (thsstp == 0.0) ? 1 : 1 + int((thslst - ths) / thsstp);
    status = hdf_file->read("PH_START", "FLOAT64_(1)", &ph);
    status = hdf_file->read("PH_STEP", "FLOAT64_(1)", &phstp);
    status = hdf_file->read("PH_END", "FLOAT64_(1)", &phlst);
    _num_phi = (phstp == 0.0) ? 1 : 1 + int((phlst - ph) / phstp);
    status = hdf_file->read("PHS_START", "FLOAT64_(1)", &phs);
    status = hdf_file->read("PHS_STEP", "FLOAT64_(1)", &phsstp);
    status = hdf_file->read("PHS_END", "FLOAT64_(1)", &phslst);
    _num_phis = (phsstp == 0.0) ? 1 : 1 + int((phslst - phs) / phsstp);
    ndirs = _num_theta * _num_thetas * _num_phi * _num_phis;
    status = hdf_file->read("EXRI", "FLOAT64_(1)", &exri);
    status = hdf_file->read("IDFC", "INT32_(1)", &idfc);
    status = hdf_file->read("XI1", "INT32_(1)", &_first_xi);
    status = hdf_file->read("NXI", "INT32_(1)", &xi_block_size);
    nxi = (_first_xi == 1) ? xi_block_size : xi_block_size + _first_xi;
    str_type = "INT32_(" + to_string(xi_block_size) + ")";
    vec_jxi = new int[xi_block_size];
    status = hdf_file->read("VEC_JXI", str_type, vec_jxi);
    str_type = "INT16_(" + to_string(xi_block_size) + ")";
    vec_ier = new short[xi_block_size];
    status = hdf_file->read("VEC_IER", str_type, vec_ier);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    vec_vk = new double[xi_block_size];
    status = hdf_file->read("VEC_VK", str_type, vec_vk);
    vec_xi = new double[xi_block_size];
    status = hdf_file->read("VEC_VK", str_type, vec_xi);
    status = hdf_file->read("NCONF", "INT32_(1)", &configurations);
    str_type = "FLOAT64_(" + to_string(xi_block_size * (configurations + 1)) + ")";
    vec_sphere_sizes = new double[xi_block_size * (configurations + 1)];
    status = hdf_file->read("VEC_SPH_SIZES", str_type, vec_sphere_sizes);
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size * (configurations + 1)) + ")";
    vec_sphere_ref_indices = new dcomplex[xi_block_size * (configurations + 1)];
    status = hdf_file->read("VEC_SPH_REFRI", str_type, vec_sphere_ref_indices);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    vec_scs1 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCS1", str_type, vec_scs1);
    vec_scs2 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCS2", str_type, vec_scs2);
    vec_abs1 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABS1", str_type, vec_abs1);
    vec_abs2 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABS2", str_type, vec_abs2);
    vec_exs1 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXS1", str_type, vec_exs1);
    vec_exs2 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXS2", str_type, vec_exs2);
    vec_albeds1 = new double[xi_block_size];
    status = hdf_file->read("VEC_ALBEDS1", str_type, vec_albeds1);
    vec_albeds2 = new double[xi_block_size];
    status = hdf_file->read("VEC_ALBEDS2", str_type, vec_albeds2);
    vec_scsrt1 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCSRT1", str_type, vec_scsrt1);
    vec_scsrt2 = new double[xi_block_size];
    status = hdf_file->read("VEC_SCSRT2", str_type, vec_scsrt2);
    vec_absrt1 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABSRT1", str_type, vec_absrt1);
    vec_absrt2 = new double[xi_block_size];
    status = hdf_file->read("VEC_ABSRT2", str_type, vec_absrt2);
    vec_exsrt1 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXSRT1", str_type, vec_exsrt1);
    vec_exsrt2 = new double[xi_block_size];
    status = hdf_file->read("VEC_EXSRT2", str_type, vec_exsrt2);
    vec_qschu1 = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCHU1", str_type, vec_qschu1);
    vec_qschu2 = new double[xi_block_size];
    status = hdf_file->read("VEC_QSCHU2", str_type, vec_qschu2);
    vec_pschu1 = new double[xi_block_size];
    status = hdf_file->read("VEC_PSCHU1", str_type, vec_pschu1);
    vec_pschu2 = new double[xi_block_size];
    status = hdf_file->read("VEC_PSCHU2", str_type, vec_pschu2);
    vec_s0mag1 = new double[xi_block_size];
    status = hdf_file->read("VEC_S0MAG1", str_type, vec_s0mag1);
    vec_s0mag2 = new double[xi_block_size];
    status = hdf_file->read("VEC_S0MAG2", str_type, vec_s0mag2);
    vec_cosav1 = new double[xi_block_size];
    status = hdf_file->read("VEC_COSAV1", str_type, vec_cosav1);
    vec_cosav2 = new double[xi_block_size];
    status = hdf_file->read("VEC_COSAV2", str_type, vec_cosav2);
    vec_raprs1 = new double[xi_block_size];
    status = hdf_file->read("VEC_RAPRS1", str_type, vec_raprs1);
    vec_raprs2 = new double[xi_block_size];
    status = hdf_file->read("VEC_RAPRS2", str_type, vec_raprs2);
    vec_fk1 = new double[xi_block_size];
    status = hdf_file->read("VEC_FK1", str_type, vec_fk1);
    vec_fk2 = new double[xi_block_size];
    status = hdf_file->read("VEC_FK2", str_type, vec_fk2);
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size) + ")";
    vec_fsas11 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAS11", str_type, vec_fsas11);
    vec_fsas21 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAS21", str_type, vec_fsas21);
    vec_fsas22 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAS22", str_type, vec_fsas22);
    vec_fsas12 = new dcomplex[xi_block_size];
    status = hdf_file->read("VEC_FSAS12", str_type, vec_fsas12);
    // Initialize directions (they are scale-independent)
    vec_dir_tidg = new double[_num_theta];
    vec_dir_tsdg = new double[_num_thetas];
    vec_dir_pidg = new double[_num_phi];
    vec_dir_psdg = new double[_num_phis];
    double cti = th, cpi = ph, cts = ths, cps = phs;
    for (int di = 0; di < _num_theta; di++) {
      vec_dir_tidg[di] = cti;
      cti += thstp;
    }
    for (int di = 0; di < _num_thetas; di++) {
      vec_dir_tsdg[di] = cts;
      cts += thsstp;
    }
    for (int di = 0; di < _num_phi; di++) {
      vec_dir_pidg[di] = cpi;
      cpi += phstp;
    }
    for (int di = 0; di < _num_phis; di++) {
      vec_dir_psdg[di] = cps;
      cps += phsstp;
    }
    str_type = "FLOAT64_(" + to_string(ndirs) + ")";
    vec_dir_scand = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SCAND", str_type, vec_dir_scand);
    vec_dir_cfmp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_CFMP", str_type, vec_dir_cfmp);
    vec_dir_sfmp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SFMP", str_type, vec_dir_sfmp);
    vec_dir_cfsp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_CFSP", str_type, vec_dir_cfsp);
    vec_dir_sfsp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SFSP", str_type, vec_dir_sfsp);
    str_type = "FLOAT64_(" + to_string(3 * ndirs) + ")";
    vec_dir_un = new double[3 * ndirs];
    status = hdf_file->read("VEC_DIR_UN", str_type, vec_dir_un);
    vec_dir_uns = new double[3 * ndirs];
    status = hdf_file->read("VEC_DIR_UNS", str_type, vec_dir_uns);
    str_type = "FLOAT64_(" + to_string(ndirs * xi_block_size) + ")";
    vec_dir_scs1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCS1", str_type, vec_dir_scs1);
    vec_dir_scs2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCS2", str_type, vec_dir_scs2);
    vec_dir_abs1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABS1", str_type, vec_dir_abs1);
    vec_dir_abs2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABS2", str_type, vec_dir_abs2);
    vec_dir_exs1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXS1", str_type, vec_dir_exs1);
    vec_dir_exs2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXS2", str_type, vec_dir_exs2);
    vec_dir_albeds1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ALBEDS1", str_type, vec_dir_albeds1);
    vec_dir_albeds2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ALBEDS2", str_type, vec_dir_albeds2);
    vec_dir_scsrt1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCSRT1", str_type, vec_dir_scsrt1);
    vec_dir_scsrt2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SCSRT2", str_type, vec_dir_scsrt2);
    vec_dir_absrt1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABSRT1", str_type, vec_dir_absrt1);
    vec_dir_absrt2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_ABSRT2", str_type, vec_dir_absrt2);
    vec_dir_exsrt1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXSRT1", str_type, vec_dir_exsrt1);
    vec_dir_exsrt2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_EXSRT2", str_type, vec_dir_exsrt2);
    str_type = "FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")";
    vec_dir_fsas11 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAS11", str_type, vec_dir_fsas11);
    vec_dir_fsas21 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAS21", str_type, vec_dir_fsas21);
    vec_dir_fsas12 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAS12", str_type, vec_dir_fsas12);
    vec_dir_fsas22 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FSAS22", str_type, vec_dir_fsas22);
    vec_dir_sas11 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS11", str_type, vec_dir_sas11);
    vec_dir_sas21 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS21", str_type, vec_dir_sas21);
    vec_dir_sas12 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS12", str_type, vec_dir_sas12);
    vec_dir_sas22 = new dcomplex[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS22", str_type, vec_dir_sas22);
    str_type = "FLOAT64_(" + to_string(ndirs * xi_block_size) + ")";
    vec_dir_qschu1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QSCHU1", str_type, vec_dir_qschu1);
    vec_dir_qschu2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_QSCHU2", str_type, vec_dir_qschu2);
    vec_dir_pschu1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_PSCHU1", str_type, vec_dir_pschu1);
    vec_dir_pschu2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_PSCHU2", str_type, vec_dir_pschu2);
    vec_dir_s0mag1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_S0MAG1", str_type, vec_dir_s0mag1);
    vec_dir_s0mag2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_S0MAG2", str_type, vec_dir_s0mag2);
    vec_dir_cosav1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_COSAV1", str_type, vec_dir_cosav1);
    vec_dir_cosav2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_COSAV2", str_type, vec_dir_cosav2);
    vec_dir_rapr1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_RAPR1", str_type, vec_dir_rapr1);
    vec_dir_rapr2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_RAPR2", str_type, vec_dir_rapr2);
    vec_dir_fl1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FL1", str_type, vec_dir_fl1);
    vec_dir_fl2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FL2", str_type, vec_dir_fl2);
    vec_dir_fr1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FR1", str_type, vec_dir_fr1);
    vec_dir_fr2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FR2", str_type, vec_dir_fr2);
    vec_dir_fk1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FK1", str_type, vec_dir_fk1);
    vec_dir_fk2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FK2", str_type, vec_dir_fk2);
    vec_dir_fx1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FX1", str_type, vec_dir_fx1);
    vec_dir_fx2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FX2", str_type, vec_dir_fx2);
    vec_dir_fy1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FY1", str_type, vec_dir_fy1);
    vec_dir_fy2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FY2", str_type, vec_dir_fy2);
    vec_dir_fz1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FZ1", str_type, vec_dir_fz1);
    vec_dir_fz2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_FZ2", str_type, vec_dir_fz2);
    vec_dir_tqel1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEL1", str_type, vec_dir_tqel1);
    vec_dir_tqel2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEL2", str_type, vec_dir_tqel2);
    vec_dir_tqer1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQER1", str_type, vec_dir_tqer1);
    vec_dir_tqer2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQER2", str_type, vec_dir_tqer2);
    vec_dir_tqek1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEK1", str_type, vec_dir_tqek1);
    vec_dir_tqek2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEK2", str_type, vec_dir_tqek2);
    vec_dir_tqex1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEX1", str_type, vec_dir_tqex1);
    vec_dir_tqex2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEX2", str_type, vec_dir_tqex2);
    vec_dir_tqey1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEY1", str_type, vec_dir_tqey1);
    vec_dir_tqey2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEY2", str_type, vec_dir_tqey2);
    vec_dir_tqez1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEZ1", str_type, vec_dir_tqez1);
    vec_dir_tqez2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQEZ2", str_type, vec_dir_tqez2);
    vec_dir_tqsl1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSL1", str_type, vec_dir_tqsl1);
    vec_dir_tqsl2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSL2", str_type, vec_dir_tqsl2);
    vec_dir_tqsr1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSR1", str_type, vec_dir_tqsr1);
    vec_dir_tqsr2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSR2", str_type, vec_dir_tqsr2);
    vec_dir_tqsk1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSK1", str_type, vec_dir_tqsk1);
    vec_dir_tqsk2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSK2", str_type, vec_dir_tqsk2);
    vec_dir_tqsx1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSX1", str_type, vec_dir_tqsx1);
    vec_dir_tqsx2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSX2", str_type, vec_dir_tqsx2);
    vec_dir_tqsy1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSY1", str_type, vec_dir_tqsy1);
    vec_dir_tqsy2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSY2", str_type, vec_dir_tqsy2);
    vec_dir_tqsz1 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSZ1", str_type, vec_dir_tqsz1);
    vec_dir_tqsz2 = new double[ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_TQSZ2", str_type, vec_dir_tqsz2);
    str_type = "FLOAT64_(" + to_string(16 * ndirs * xi_block_size) + ")";
    vec_dir_mull = new double[16 * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULL", str_type, vec_dir_mull);
    vec_dir_mulllr = new double[16 * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULLLR", str_type, vec_dir_mulllr);
    status = hdf_file->close();
    delete hdf_file;
  } else {
    if (hdf_file != NULL) delete hdf_file;
    UnrecognizedFormatException ex("Error: " + hdf5_name + " not recognized as a valid HDF5 file!");
    throw ex;
  }
}

InclusionOutputInfo::InclusionOutputInfo(const int flag) {
  /*
    create a dummy placeholder just to know I should skip MPI_Send and MPI_Recv
  */
  if (flag == 1) {
    _skip_flag = 1;
  } else {
    UnrecognizedOutputInfo ex(flag);
    throw ex;
  }
}

InclusionOutputInfo::~InclusionOutputInfo() {
  if (_skip_flag != 1) {
    delete[] vec_x_coords;
    delete[] vec_y_coords;
    delete[] vec_z_coords;
    delete[] vec_jxi;
    delete[] vec_ier;
    delete[] vec_vk;
    delete[] vec_xi;
    delete[] vec_sphere_sizes;
    delete[] vec_sphere_ref_indices;
    delete[] vec_scs1;
    delete[] vec_scs2;
    delete[] vec_abs1;
    delete[] vec_abs2;
    delete[] vec_exs1;
    delete[] vec_exs2;
    delete[] vec_albeds1;
    delete[] vec_albeds2;
    delete[] vec_scsrt1;
    delete[] vec_scsrt2;
    delete[] vec_absrt1;
    delete[] vec_absrt2;
    delete[] vec_exsrt1;
    delete[] vec_exsrt2;
    delete[] vec_qschu1;
    delete[] vec_qschu2;
    delete[] vec_pschu1;
    delete[] vec_pschu2;
    delete[] vec_s0mag1;
    delete[] vec_s0mag2;
    delete[] vec_cosav1;
    delete[] vec_cosav2;
    delete[] vec_raprs1;
    delete[] vec_raprs2;
    delete[] vec_fk1;
    delete[] vec_fk2;
    delete[] vec_fsas11;
    delete[] vec_fsas21;
    delete[] vec_fsas22;
    delete[] vec_fsas12;
    delete[] vec_dir_tidg;
    delete[] vec_dir_pidg;
    delete[] vec_dir_tsdg;
    delete[] vec_dir_psdg;
    delete[] vec_dir_scand;
    delete[] vec_dir_cfmp;
    delete[] vec_dir_sfmp;
    delete[] vec_dir_cfsp;
    delete[] vec_dir_sfsp;
    delete[] vec_dir_un;
    delete[] vec_dir_uns;
    delete[] vec_dir_scs1;
    delete[] vec_dir_scs2;
    delete[] vec_dir_abs1;
    delete[] vec_dir_abs2;
    delete[] vec_dir_exs1;
    delete[] vec_dir_exs2;
    delete[] vec_dir_albeds1;
    delete[] vec_dir_albeds2;
    delete[] vec_dir_scsrt1;
    delete[] vec_dir_scsrt2;
    delete[] vec_dir_absrt1;
    delete[] vec_dir_absrt2;
    delete[] vec_dir_exsrt1;
    delete[] vec_dir_exsrt2;
    delete[] vec_dir_fsas11;
    delete[] vec_dir_fsas21;
    delete[] vec_dir_fsas12;
    delete[] vec_dir_fsas22;
    delete[] vec_dir_sas11;
    delete[] vec_dir_sas21;
    delete[] vec_dir_sas12;
    delete[] vec_dir_sas22;
    delete[] vec_dir_qschu1;
    delete[] vec_dir_qschu2;
    delete[] vec_dir_pschu1;
    delete[] vec_dir_pschu2;
    delete[] vec_dir_s0mag1;
    delete[] vec_dir_s0mag2;
    delete[] vec_dir_cosav1;
    delete[] vec_dir_cosav2;
    delete[] vec_dir_rapr1;
    delete[] vec_dir_rapr2;
    delete[] vec_dir_fl1;
    delete[] vec_dir_fl2;
    delete[] vec_dir_fr1;
    delete[] vec_dir_fr2;
    delete[] vec_dir_fk1;
    delete[] vec_dir_fk2;
    delete[] vec_dir_fx1;
    delete[] vec_dir_fx2;
    delete[] vec_dir_fy1;
    delete[] vec_dir_fy2;
    delete[] vec_dir_fz1;
    delete[] vec_dir_fz2;
    delete[] vec_dir_tqel1;
    delete[] vec_dir_tqel2;
    delete[] vec_dir_tqer1;
    delete[] vec_dir_tqer2;
    delete[] vec_dir_tqek1;
    delete[] vec_dir_tqek2;
    delete[] vec_dir_tqex1;
    delete[] vec_dir_tqex2;
    delete[] vec_dir_tqey1;
    delete[] vec_dir_tqey2;
    delete[] vec_dir_tqez1;
    delete[] vec_dir_tqez2;
    delete[] vec_dir_tqsl1;
    delete[] vec_dir_tqsl2;
    delete[] vec_dir_tqsr1;
    delete[] vec_dir_tqsr2;
    delete[] vec_dir_tqsk1;
    delete[] vec_dir_tqsk2;
    delete[] vec_dir_tqsx1;
    delete[] vec_dir_tqsx2;
    delete[] vec_dir_tqsy1;
    delete[] vec_dir_tqsy2;
    delete[] vec_dir_tqsz1;
    delete[] vec_dir_tqsz2;
    delete[] vec_dir_mull;
    delete[] vec_dir_mulllr;
  }
}

long InclusionOutputInfo::compute_size() {
  long result = sizeof(np_int);
  result += 21 * sizeof(int);
  result += 13 * sizeof(double);
  result += 120 * sizeof(long); // vector pointers
  result += xi_block_size * (sizeof(int) + sizeof(short)); // vec_jxi, vec_ier
  result += 3 * nsph * sizeof(double); // coordinate vectors
  result += 29 * xi_block_size * sizeof(double); // wavelength dependent real values
  result += 5 * xi_block_size * sizeof(dcomplex); // wavelength dependent complex values
  result += 15 * ndirs * sizeof(double); // values depending only on directions
  result += 92 * ndirs * xi_block_size * sizeof(double); // values depending on directions and wavelengths
  result += 8 * ndirs * xi_block_size * sizeof(dcomplex); // values depending on directions and wavelengths
  return result;
}

long InclusionOutputInfo::compute_size(
  ScattererConfiguration *sc, GeometryConfiguration *gc,
  int first_xi, int xi_length
) {
  long result = sizeof(np_int);
  result += 21 * sizeof(int);
  result += 13 * sizeof(double);
  result += 120 * sizeof(long); // vector pointers
  int _nsph = gc->number_of_spheres;
  double _th = gc->in_theta_start;
  double _thstp = gc->in_theta_step;
  double _thlst = gc->in_theta_end;
  double _ths = gc->sc_theta_start;
  double _thsstp = gc->sc_theta_step;
  double _thslst = gc->sc_theta_end;
  int num_theta = 1 + int((_thlst - _th) / _thstp);
  int num_thetas = 1 + int((_thslst - _ths) / _thsstp);
  double _ph = gc->in_phi_start;
  double _phstp = gc->in_phi_step;
  double _phlst = gc->in_phi_end;
  double _phs = gc->sc_phi_start;
  double _phsstp = gc->sc_phi_step;
  double _phslst = gc->sc_phi_end;
  int num_phi = 1 + int((_phlst - _ph) / _phstp);
  int num_phis = 1 + int((_phslst - _phs) / _phsstp);
  int _ndirs = num_theta * num_thetas * num_phi * num_phis;
  int _nxi = sc->number_of_scales;
  int _xi_block_size = (xi_length == 0) ? _nxi : xi_length;
  int _configurations = sc->configurations;
  result += _xi_block_size * (sizeof(int) + sizeof(short)); // vec_jxi, vec_ier
  result += 3 * _nsph * sizeof(double); // coordinate vectors
  result += 29 * _xi_block_size * sizeof(double); // wavelength dependent real values
  result += 5 * _xi_block_size * sizeof(dcomplex); // wavelength dependent complex values
  result += 15 * _ndirs * sizeof(double); // values depending only on directions
  result += 92 * _ndirs * _xi_block_size * sizeof(double); // values depending on directions and wavelengths
  result += 8 * _ndirs * _xi_block_size * sizeof(dcomplex); // values depending on directions and wavelengths
  return result;
}

int InclusionOutputInfo::insert(const InclusionOutputInfo &rhs) {
  int result = 0;
  if (rhs.skip_flag != 1) {
    result += (rhs.nsph == nsph) ? 0 : 1;
    result += (rhs.inpol == inpol) ? 0 : 1;
    result += (rhs.iavm == iavm) ? 0 : 1;
    result += (rhs.isam == isam) ? 0 : 1;
    result += (rhs._num_theta == _num_theta) ? 0 : 1;
    result += (rhs._num_thetas == _num_thetas) ? 0 : 1;
    result += (rhs._num_phi == _num_phi) ? 0 : 1;
    result += (rhs._num_phis == _num_phis) ? 0 : 1;
    result += (rhs.ndirs == ndirs) ? 0 : 1;
    result += (rhs.exri == exri) ? 0 : 1;
    result += (rhs.idfc == idfc) ? 0 : 1;
    result += (rhs.configurations == configurations) ? 0 : 1;
    if (result == 0) {
      int offset, chunk_size, xi1;
      xi1 = rhs._first_xi;
      // Insert vectors whose size depends on wavelengths
      offset = xi1 - _first_xi;
      chunk_size = rhs.xi_block_size;
      memcpy(vec_jxi + offset, rhs.vec_jxi, chunk_size * sizeof(int));
      memcpy(vec_ier + offset, rhs.vec_ier, chunk_size * sizeof(short));
      memcpy(vec_vk + offset, rhs.vec_vk, chunk_size * sizeof(double));
      memcpy(vec_xi + offset, rhs.vec_xi, chunk_size * sizeof(double));
      memcpy(vec_scs1 + offset, rhs.vec_scs1, chunk_size * sizeof(double));
      memcpy(vec_scs2 + offset, rhs.vec_scs2, chunk_size * sizeof(double));
      memcpy(vec_abs1 + offset, rhs.vec_abs1, chunk_size * sizeof(double));
      memcpy(vec_abs2 + offset, rhs.vec_abs2, chunk_size * sizeof(double));
      memcpy(vec_exs1 + offset, rhs.vec_exs1, chunk_size * sizeof(double));
      memcpy(vec_exs2 + offset, rhs.vec_exs2, chunk_size * sizeof(double));
      memcpy(vec_albeds1 + offset, rhs.vec_albeds1, chunk_size * sizeof(double));
      memcpy(vec_albeds2 + offset, rhs.vec_albeds2, chunk_size * sizeof(double));
      memcpy(vec_scsrt1 + offset, rhs.vec_scsrt1, chunk_size * sizeof(double));
      memcpy(vec_scsrt2 + offset, rhs.vec_scsrt2, chunk_size * sizeof(double));
      memcpy(vec_absrt1 + offset, rhs.vec_absrt1, chunk_size * sizeof(double));
      memcpy(vec_absrt2 + offset, rhs.vec_absrt2, chunk_size * sizeof(double));
      memcpy(vec_exsrt1 + offset, rhs.vec_exsrt1, chunk_size * sizeof(double));
      memcpy(vec_exsrt2 + offset, rhs.vec_exsrt2, chunk_size * sizeof(double));
      memcpy(vec_qschu1 + offset, rhs.vec_qschu1, chunk_size * sizeof(double));
      memcpy(vec_qschu2 + offset, rhs.vec_qschu2, chunk_size * sizeof(double));
      memcpy(vec_pschu1 + offset, rhs.vec_pschu1, chunk_size * sizeof(double));
      memcpy(vec_pschu2 + offset, rhs.vec_pschu2, chunk_size * sizeof(double));
      memcpy(vec_s0mag1 + offset, rhs.vec_s0mag1, chunk_size * sizeof(double));
      memcpy(vec_s0mag2 + offset, rhs.vec_s0mag2, chunk_size * sizeof(double));
      memcpy(vec_cosav1 + offset, rhs.vec_cosav1, chunk_size * sizeof(double));
      memcpy(vec_cosav2 + offset, rhs.vec_cosav2, chunk_size * sizeof(double));
      memcpy(vec_raprs1 + offset, rhs.vec_raprs1, chunk_size * sizeof(double));
      memcpy(vec_raprs2 + offset, rhs.vec_raprs2, chunk_size * sizeof(double));
      memcpy(vec_fk1 + offset, rhs.vec_fk1, chunk_size * sizeof(double));
      memcpy(vec_fk2 + offset, rhs.vec_fk2, chunk_size * sizeof(double));
      memcpy(vec_fsas11 + offset, rhs.vec_fsas11, chunk_size * sizeof(dcomplex));
      memcpy(vec_fsas21 + offset, rhs.vec_fsas21, chunk_size * sizeof(dcomplex));
      memcpy(vec_fsas22 + offset, rhs.vec_fsas22, chunk_size * sizeof(dcomplex));
      memcpy(vec_fsas12 + offset, rhs.vec_fsas12, chunk_size * sizeof(dcomplex));
    
      // Insert vectors whose sizes depend on configurations
      offset = (xi1 - _first_xi) * (configurations + 1);
      chunk_size = rhs.xi_block_size * (configurations + 1);
      memcpy(vec_sphere_sizes + offset, rhs.vec_sphere_sizes, chunk_size * sizeof(double));
      memcpy(vec_sphere_ref_indices + offset, rhs.vec_sphere_ref_indices, chunk_size * sizeof(dcomplex));

      // Insert vectors whose sizes depend on directions and wavelengths
      offset = (xi1 - _first_xi) * ndirs;
      chunk_size = rhs.xi_block_size * ndirs;
      memcpy(vec_dir_scs1 + offset, rhs.vec_dir_scs1, chunk_size * sizeof(double));
      memcpy(vec_dir_scs2 + offset, rhs.vec_dir_scs2, chunk_size * sizeof(double));
      memcpy(vec_dir_abs1 + offset, rhs.vec_dir_abs1, chunk_size * sizeof(double));
      memcpy(vec_dir_abs2 + offset, rhs.vec_dir_abs2, chunk_size * sizeof(double));
      memcpy(vec_dir_exs1 + offset, rhs.vec_dir_exs1, chunk_size * sizeof(double));
      memcpy(vec_dir_exs2 + offset, rhs.vec_dir_exs2, chunk_size * sizeof(double));
      memcpy(vec_dir_albeds1 + offset, rhs.vec_dir_albeds1, chunk_size * sizeof(double));
      memcpy(vec_dir_albeds2 + offset, rhs.vec_dir_albeds2, chunk_size * sizeof(double));
      memcpy(vec_dir_scsrt1 + offset, rhs.vec_dir_scsrt1, chunk_size * sizeof(double));
      memcpy(vec_dir_scsrt2 + offset, rhs.vec_dir_scsrt2, chunk_size * sizeof(double));
      memcpy(vec_dir_absrt1 + offset, rhs.vec_dir_absrt1, chunk_size * sizeof(double));
      memcpy(vec_dir_absrt2 + offset, rhs.vec_dir_absrt2, chunk_size * sizeof(double));
      memcpy(vec_dir_exsrt1 + offset, rhs.vec_dir_exsrt1, chunk_size * sizeof(double));
      memcpy(vec_dir_exsrt2 + offset, rhs.vec_dir_exsrt2, chunk_size * sizeof(double));
      memcpy(vec_dir_fsas11 + offset, rhs.vec_dir_fsas11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_fsas21 + offset, rhs.vec_dir_fsas21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_fsas12 + offset, rhs.vec_dir_fsas12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_fsas22 + offset, rhs.vec_dir_fsas22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas11 + offset, rhs.vec_dir_sas11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas21 + offset, rhs.vec_dir_sas21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas12 + offset, rhs.vec_dir_sas12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas22 + offset, rhs.vec_dir_sas22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_qschu1 + offset, rhs.vec_dir_qschu1, chunk_size * sizeof(double));
      memcpy(vec_dir_qschu2 + offset, rhs.vec_dir_qschu2, chunk_size * sizeof(double));
      memcpy(vec_dir_pschu1 + offset, rhs.vec_dir_pschu1, chunk_size * sizeof(double));
      memcpy(vec_dir_pschu2 + offset, rhs.vec_dir_pschu2, chunk_size * sizeof(double));
      memcpy(vec_dir_s0mag1 + offset, rhs.vec_dir_s0mag1, chunk_size * sizeof(double));
      memcpy(vec_dir_s0mag2 + offset, rhs.vec_dir_s0mag2, chunk_size * sizeof(double));
      memcpy(vec_dir_cosav1 + offset, rhs.vec_dir_cosav1, chunk_size * sizeof(double));
      memcpy(vec_dir_cosav2 + offset, rhs.vec_dir_cosav2, chunk_size * sizeof(double));
      memcpy(vec_dir_rapr1 + offset, rhs.vec_dir_rapr1, chunk_size * sizeof(double));
      memcpy(vec_dir_rapr2 + offset, rhs.vec_dir_rapr2, chunk_size * sizeof(double));
      memcpy(vec_dir_fl1 + offset, rhs.vec_dir_fl1, chunk_size * sizeof(double));
      memcpy(vec_dir_fl2 + offset, rhs.vec_dir_fl2, chunk_size * sizeof(double));
      memcpy(vec_dir_fr1 + offset, rhs.vec_dir_fr1, chunk_size * sizeof(double));
      memcpy(vec_dir_fr2 + offset, rhs.vec_dir_fr2, chunk_size * sizeof(double));
      memcpy(vec_dir_fk1 + offset, rhs.vec_dir_fk1, chunk_size * sizeof(double));
      memcpy(vec_dir_fk2 + offset, rhs.vec_dir_fk2, chunk_size * sizeof(double));
      memcpy(vec_dir_fx1 + offset, rhs.vec_dir_fx1, chunk_size * sizeof(double));
      memcpy(vec_dir_fx2 + offset, rhs.vec_dir_fx2, chunk_size * sizeof(double));
      memcpy(vec_dir_fy1 + offset, rhs.vec_dir_fy1, chunk_size * sizeof(double));
      memcpy(vec_dir_fy2 + offset, rhs.vec_dir_fy2, chunk_size * sizeof(double));
      memcpy(vec_dir_fz1 + offset, rhs.vec_dir_fz1, chunk_size * sizeof(double));
      memcpy(vec_dir_fz2 + offset, rhs.vec_dir_fz2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqel1 + offset, rhs.vec_dir_tqel1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqel2 + offset, rhs.vec_dir_tqel2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqer1 + offset, rhs.vec_dir_tqer1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqer2 + offset, rhs.vec_dir_tqer2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqek1 + offset, rhs.vec_dir_tqek1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqek2 + offset, rhs.vec_dir_tqek2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqex1 + offset, rhs.vec_dir_tqex1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqex2 + offset, rhs.vec_dir_tqex2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqey1 + offset, rhs.vec_dir_tqey1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqey2 + offset, rhs.vec_dir_tqey2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqez1 + offset, rhs.vec_dir_tqez1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqez2 + offset, rhs.vec_dir_tqez2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsl1 + offset, rhs.vec_dir_tqsl1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsl2 + offset, rhs.vec_dir_tqsl2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsr1 + offset, rhs.vec_dir_tqsr1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsr2 + offset, rhs.vec_dir_tqsr2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsk1 + offset, rhs.vec_dir_tqsk1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsk2 + offset, rhs.vec_dir_tqsk2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsx1 + offset, rhs.vec_dir_tqsx1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsx2 + offset, rhs.vec_dir_tqsx2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsy1 + offset, rhs.vec_dir_tqsy1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsy2 + offset, rhs.vec_dir_tqsy2, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsz1 + offset, rhs.vec_dir_tqsz1, chunk_size * sizeof(double));
      memcpy(vec_dir_tqsz2 + offset, rhs.vec_dir_tqsz2, chunk_size * sizeof(double));
      memcpy(vec_dir_mull + 16 * offset, rhs.vec_dir_mull, 16 * chunk_size * sizeof(double));
      memcpy(vec_dir_mulllr + 16 * offset, rhs.vec_dir_mulllr, 16 * chunk_size * sizeof(double));
    }
  }
  return result;
}

int InclusionOutputInfo::write(const std::string &output, const std::string &format) {
  int result = 0;
  if (format.compare("LEGACY") == 0) {
    result = write_legacy(output);
  } else if (format.compare("HDF5") == 0) {
    result = write_hdf5(output);
  } else {
    string message = "Unknown format mode: \"" + format + "\"";
    throw UnrecognizedConfigurationException(message);
  }
  return result;
}

int InclusionOutputInfo::write_hdf5(const std::string &file_name) {
  List<string> *rec_name_list = new List<string>(1);
  List<string> *rec_type_list = new List<string>(1);
  List<void *> *rec_ptr_list = new List<void *>(1);
  string str_type, str_name;
  rec_name_list->set(0, "NSPH");
  rec_type_list->set(0, "INT32_(1)");
  rec_ptr_list->set(0, &nsph);
  rec_name_list->append("LI");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&li);
  rec_name_list->append("LE");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&le);
  rec_name_list->append("LM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&lm);
  rec_name_list->append("MXNDM");
  rec_type_list->append("INT64_(1)");
  rec_ptr_list->append(&mxndm);
  rec_name_list->append("INPOL");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&inpol);
  rec_name_list->append("NPNT");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&npnt);
  rec_name_list->append("NPNTTS");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&npntts);
  rec_name_list->append("IAVM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&iavm);
  rec_name_list->append("ISAM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&isam);
  rec_name_list->append("JWTM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&jwtm);
  rec_name_list->append("VEC_SPH_X");
  rec_type_list->append("FLOAT64_(" + to_string(nsph) + ")");
  rec_ptr_list->append(vec_x_coords);
  rec_name_list->append("VEC_SPH_Y");
  rec_type_list->append("FLOAT64_(" + to_string(nsph) + ")");
  rec_ptr_list->append(vec_y_coords);
  rec_name_list->append("VEC_SPH_Z");
  rec_type_list->append("FLOAT64_(" + to_string(nsph) + ")");
  rec_ptr_list->append(vec_z_coords);
  rec_name_list->append("TH_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&th);
  rec_name_list->append("TH_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thstp);
  rec_name_list->append("TH_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thlst);
  rec_name_list->append("THS_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&ths);
  rec_name_list->append("THS_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thsstp);
  rec_name_list->append("THS_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thslst);
  rec_name_list->append("PH_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&ph);
  rec_name_list->append("PH_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phstp);
  rec_name_list->append("PH_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phlst);
  rec_name_list->append("PHS_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phs);
  rec_name_list->append("PHS_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phsstp);
  rec_name_list->append("PHS_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phslst);
  rec_name_list->append("EXRI");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&exri);
  rec_name_list->append("IDFC");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&idfc);
  rec_name_list->append("XI1");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&_first_xi);
  rec_name_list->append("NXI");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&xi_block_size);
  rec_name_list->append("VEC_JXI");
  rec_type_list->append("INT32_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_jxi);
  rec_name_list->append("VEC_IER");
  rec_type_list->append("INT16_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_ier);
  rec_name_list->append("VEC_VK");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_vk);
  rec_name_list->append("VEC_XI");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_xi);
  rec_name_list->append("NCONF");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&configurations);
  rec_name_list->append("VEC_SPH_SIZES");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size * (configurations + 1)) + ")");
  rec_ptr_list->append(vec_sphere_sizes);
  rec_name_list->append("VEC_SPH_REFRI");
  rec_type_list->append("FLOAT64_(" + to_string(2 * xi_block_size * (configurations+ 1)) + ")");
  rec_ptr_list->append(vec_sphere_ref_indices);
  str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
  rec_name_list->append("VEC_SCS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_scs1);
  rec_name_list->append("VEC_SCS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_scs2);
  rec_name_list->append("VEC_ABS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_abs1);
  rec_name_list->append("VEC_ABS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_abs2);
  rec_name_list->append("VEC_EXS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_exs1);
  rec_name_list->append("VEC_EXS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_exs2);
  rec_name_list->append("VEC_ALBEDS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_albeds1);
  rec_name_list->append("VEC_ALBEDS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_albeds2);
  rec_name_list->append("VEC_SCSRT1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_scsrt1);
  rec_name_list->append("VEC_SCSRT2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_scsrt2);
  rec_name_list->append("VEC_ABSRT1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_absrt1);
  rec_name_list->append("VEC_ABSRT2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_absrt2);
  rec_name_list->append("VEC_EXSRT1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_exsrt1);
  rec_name_list->append("VEC_EXSRT2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_exsrt2);
  rec_name_list->append("VEC_QSCHU1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_qschu1);
  rec_name_list->append("VEC_QSCHU2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_qschu2);
  rec_name_list->append("VEC_PSCHU1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_pschu1);
  rec_name_list->append("VEC_PSCHU2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_pschu2);
  rec_name_list->append("VEC_S0MAG1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_s0mag1);
  rec_name_list->append("VEC_S0MAG2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_s0mag2);
  rec_name_list->append("VEC_COSAV1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_cosav1);
  rec_name_list->append("VEC_COSAV2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_cosav2);
  rec_name_list->append("VEC_RAPRS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_raprs1);
  rec_name_list->append("VEC_RAPRS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_raprs2);
  rec_name_list->append("VEC_FK1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fk1);
  rec_name_list->append("VEC_FK2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fk2);
  str_type = "FLOAT64_(" + to_string(2 * xi_block_size) + ")";
  rec_name_list->append("VEC_FSAS11");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fsas11);
  rec_name_list->append("VEC_FSAS21");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fsas21);
  rec_name_list->append("VEC_FSAS12");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fsas12);
  rec_name_list->append("VEC_FSAS22");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fsas22);
  str_type = "FLOAT64_(" + to_string(ndirs) + ")";
  rec_name_list->append("VEC_DIR_SCAND");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_scand);
  rec_name_list->append("VEC_DIR_CFMP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_cfmp);
  rec_name_list->append("VEC_DIR_CFSP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_cfsp);
  rec_name_list->append("VEC_DIR_SFMP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sfmp);
  rec_name_list->append("VEC_DIR_SFSP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sfsp);
  str_type = "FLOAT64_(" + to_string(3 * ndirs) + ")";
  rec_name_list->append("VEC_DIR_UN");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_un);
  rec_name_list->append("VEC_DIR_UNS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_uns);
  str_type = "FLOAT64_(" + to_string(ndirs * xi_block_size) + ")";
  rec_name_list->append("VEC_DIR_SCS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_scs1);
  rec_name_list->append("VEC_DIR_SCS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_scs2);
  rec_name_list->append("VEC_DIR_ABS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_abs1);
  rec_name_list->append("VEC_DIR_ABS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_abs2);
  rec_name_list->append("VEC_DIR_EXS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_exs1);
  rec_name_list->append("VEC_DIR_EXS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_exs2);
  rec_name_list->append("VEC_DIR_ALBEDS1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_albeds1);
  rec_name_list->append("VEC_DIR_ALBEDS2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_albeds2);
  rec_name_list->append("VEC_DIR_SCSRT1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_scsrt1);
  rec_name_list->append("VEC_DIR_SCSRT2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_scsrt2);
  rec_name_list->append("VEC_DIR_ABSRT1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_absrt1);
  rec_name_list->append("VEC_DIR_ABSRT2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_absrt2);
  rec_name_list->append("VEC_DIR_EXSRT1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_exsrt1);
  rec_name_list->append("VEC_DIR_EXSRT2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_exsrt2);
  str_type = "FLOAT64_(" + to_string(2 * ndirs * xi_block_size) + ")"; 
  rec_name_list->append("VEC_DIR_FSAS11");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fsas11);
  rec_name_list->append("VEC_DIR_FSAS21");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fsas21);
  rec_name_list->append("VEC_DIR_FSAS12");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fsas12);
  rec_name_list->append("VEC_DIR_FSAS22");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fsas22);
  rec_name_list->append("VEC_DIR_SAS11");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas11);
  rec_name_list->append("VEC_DIR_SAS21");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas21);
  rec_name_list->append("VEC_DIR_SAS12");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas12);
  rec_name_list->append("VEC_DIR_SAS22");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas22);
  str_type = "FLOAT64_(" + to_string(ndirs * xi_block_size) + ")";
  rec_name_list->append("VEC_DIR_QSCHU1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_qschu1);
  rec_name_list->append("VEC_DIR_QSCHU2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_qschu2);
  rec_name_list->append("VEC_DIR_PSCHU1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_pschu1);
  rec_name_list->append("VEC_DIR_PSCHU2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_pschu2);
  rec_name_list->append("VEC_DIR_S0MAG1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_s0mag1);
  rec_name_list->append("VEC_DIR_S0MAG2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_s0mag2);
  rec_name_list->append("VEC_DIR_COSAV1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_cosav1);
  rec_name_list->append("VEC_DIR_COSAV2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_cosav2);
  rec_name_list->append("VEC_DIR_RAPR1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_rapr1);
  rec_name_list->append("VEC_DIR_RAPR2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_rapr2);
  rec_name_list->append("VEC_DIR_FL1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fl1);
  rec_name_list->append("VEC_DIR_FL2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fl2);
  rec_name_list->append("VEC_DIR_FR1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fr1);
  rec_name_list->append("VEC_DIR_FR2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fr2);
  rec_name_list->append("VEC_DIR_FK1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fk1);
  rec_name_list->append("VEC_DIR_FK2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fk2);
  rec_name_list->append("VEC_DIR_FX1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fx1);
  rec_name_list->append("VEC_DIR_FX2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fx2);
  rec_name_list->append("VEC_DIR_FY1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fy1);
  rec_name_list->append("VEC_DIR_FY2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fy2);
  rec_name_list->append("VEC_DIR_FZ1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fz1);
  rec_name_list->append("VEC_DIR_FZ2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fz2);
  rec_name_list->append("VEC_DIR_TQEL1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqel1);
  rec_name_list->append("VEC_DIR_TQEL2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqel2);
  rec_name_list->append("VEC_DIR_TQER1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqer1);
  rec_name_list->append("VEC_DIR_TQER2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqer2);
  rec_name_list->append("VEC_DIR_TQEK1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqek1);
  rec_name_list->append("VEC_DIR_TQEK2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqek2);
  rec_name_list->append("VEC_DIR_TQEX1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqex1);
  rec_name_list->append("VEC_DIR_TQEX2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqex2);
  rec_name_list->append("VEC_DIR_TQEY1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqey1);
  rec_name_list->append("VEC_DIR_TQEY2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqey2);
  rec_name_list->append("VEC_DIR_TQEZ1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqez1);
  rec_name_list->append("VEC_DIR_TQEZ2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqez2);
  rec_name_list->append("VEC_DIR_TQSL1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsl1);
  rec_name_list->append("VEC_DIR_TQSL2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsl2);
  rec_name_list->append("VEC_DIR_TQSR1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsr1);
  rec_name_list->append("VEC_DIR_TQSR2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsr2);
  rec_name_list->append("VEC_DIR_TQSK1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsk1);
  rec_name_list->append("VEC_DIR_TQSK2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsk2);
  rec_name_list->append("VEC_DIR_TQSX1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsx1);
  rec_name_list->append("VEC_DIR_TQSX2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsx2);
  rec_name_list->append("VEC_DIR_TQSY1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsy1);
  rec_name_list->append("VEC_DIR_TQSY2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsy2);
  rec_name_list->append("VEC_DIR_TQSZ1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsz1);
  rec_name_list->append("VEC_DIR_TQSZ2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_tqsz2);
  str_type = "FLOAT64_(" + to_string(16 * ndirs * xi_block_size) + ")";
  rec_name_list->append("VEC_DIR_MULL");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_mull);
  rec_name_list->append("VEC_DIR_MULLLR");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_mulllr);
  
  // Convert the lists to arrays and write them to HDF5
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
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete schema;
  delete hdf_file;
  return 0;
}

int InclusionOutputInfo::write_legacy(const std::string &output) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  int result = 0;
  FILE *p_outfile = fopen(output.c_str(), "w");
  if (p_outfile != NULL) {
    if (vec_jxi[0] == 1) {
      fprintf(p_outfile, " READ(IR,*)NSPH,LI,LE,MXNDM,INPOL,NPNT,NPNTTS,IAVM,ISAM\n");
#ifdef USE_ILP64
      fprintf(
	p_outfile, " %5d%5d%5d%5ld%5d%5d%5d%5d%5d\n",
	nsph, li, le, mxndm, inpol, npnt, npntts,
	iavm, iavm
      );
#else
      fprintf(
	p_outfile, " %5d%5d%5d%5d%5d%5d%5d%5d%5d\n",
	nsph, li, le, mxndm, inpol, npnt, npntts,
	iavm, iavm
      );
#endif // USE_ILP64
      fprintf(p_outfile, " READ(IR,*)RXX(I),RYY(I),RZZ(I)\n");
      for (int ri = 0; ri < nsph; ri++) {
	fprintf(
	  p_outfile, "%17.8lE%17.8lE%17.8lE\n",
	  vec_x_coords[ri], vec_y_coords[ri], vec_z_coords[ri]
	);
      }
      fprintf(p_outfile, " READ(IR,*)TH,THSTP,THLST,THS,THSSTP,THSLST\n");
      fprintf(
        p_outfile, " %10.3lE%10.3lE%10.3lE%10.3lE%10.3lE%10.3lE\n",
	th, thstp, thlst, ths, thsstp, thslst
      );
      fprintf(p_outfile, " READ(IR,*)PH,PHSTP,PHLST,PHS,PHSSTP,PHSLST\n");
      fprintf(
        p_outfile, " %10.3lE%10.3lE%10.3lE%10.3lE%10.3lE%10.3lE\n",
	ph, phstp, phlst, phs, phsstp, phslst
      );
      fprintf(p_outfile, " READ(IR,*)JWTM\n");
      fprintf(p_outfile, " %5d\n", jwtm);
      fprintf(p_outfile, "  READ(ITIN)NSPHT\n");
      fprintf(p_outfile, "  READ(ITIN)(IOG(I),I=1,NSPH)\n");
      fprintf(p_outfile, "  READ(ITIN)EXDC,WP,XIP,IDFC,NXI\n");
      fprintf(p_outfile, "  READ(ITIN)(XIV(I),I=1,NXI)\n");
      fprintf(p_outfile, "  READ(ITIN)NSHL(I),ROS(I)\n");
      fprintf(p_outfile, "  READ(ITIN)(RCF(I,NS),NS=1,NSH)\n \n");
      fprintf(p_outfile, "  REFR. INDEX OF EXTERNAL MEDIUM=%15.7lE\n", exri);
      if (idfc < 0) {
	fprintf( p_outfile, "  VK=%15.7lE, XI IS SCALE FACTOR FOR LENGTHS\n \n", vec_vk[0]);
      }
      // End preamble writing
    }
    // Wavelength loop
    for (int jxi = 0; jxi < xi_block_size; jxi++) {
	int done_dirs = 0;
	double alamb = 2.0 * 3.141592653589793238 / vec_vk[jxi];
	fprintf(p_outfile, "========== JXI =%3d ====================\n", vec_jxi[jxi]);
	if (idfc >= 0) {
	  fprintf(p_outfile, "  VK=%15.7lE, XI=%15.7lE\n", vec_vk[jxi], vec_xi[jxi]);
	} else {
	  fprintf(p_outfile, "  XI=%15.7lE\n", vec_xi[jxi]);
	}
	if (vec_ier[jxi] == 1) {
	  fprintf(p_outfile, "  STOP IN INDME\n");
	  break;
	} else if (vec_ier[jxi] == 2) {
	  fprintf(p_outfile, "  STOP IN OSPV\n");
	  break;
	}
	for (int i168 = 1; i168 <= configurations; i168++) {
	  int cindex = jxi * (configurations + 1) + i168 - 1;
	  if (vec_sphere_ref_indices[cindex] == cc0) {
	    fprintf(p_outfile, "  SPHERE N.%2d: SIZE=%15.7lE\n", i168, vec_sphere_sizes[cindex]);
	  } else {
	    fprintf(
	      p_outfile, "  SPHERE N.%2d: SIZE=%15.7lE, REFRACTIVE INDEX=%15.7lE%15.7lE\n",
	      i168, vec_sphere_sizes[cindex], real(vec_sphere_ref_indices[cindex]),
	      imag(vec_sphere_ref_indices[cindex])
	    );
	  }
	} // i168 configuration loop
	int cindex = jxi * (configurations + 1) + configurations;
	fprintf(
	  p_outfile, "  EXT. SPHERE: SIZE=%15.7lE, REFRACTIVE INDEX=%15.7lE%15.7lE\n",
	  vec_sphere_sizes[cindex], real(vec_sphere_ref_indices[cindex]),
	  imag(vec_sphere_ref_indices[cindex])
	);
	fprintf(p_outfile, "     ENSEMBLE AVERAGE, MODE%2d\n", iavm);
	if (inpol == 0) fprintf(p_outfile, "   LIN -1\n");
	else fprintf(p_outfile, "  CIRC -1\n");
	fprintf(p_outfile, " ----- SCS ----- ABS ----- EXS ----- ALBEDS --\n");
	fprintf(
	  p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
	  vec_scs1[jxi], vec_abs1[jxi], vec_exs1[jxi], vec_albeds1[jxi]
	);
	// fprintf(p_outfile, "INSERTION: SCASECM %5d%15.7E%15.7E%15.7E%15.7E\n",
	// 	-1, alamb, vec_scs1[jxi], vec_abs1[jxi], vec_exs1[jxi]
	// );
	fprintf(p_outfile, " ---- SCS/GS -- ABC/GS -- EXS/GS ---\n");
	fprintf(
	  p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	  vec_scsrt1[jxi], vec_absrt1[jxi], vec_exsrt1[jxi]
	);
	fprintf(
	  p_outfile, "  FSAS(1,1)=%15.7lE%15.7lE   FSAS(2,1)=%15.7lE%15.7lE\n",
	  real(vec_fsas11[jxi]), imag(vec_fsas11[jxi]),
	  real(vec_fsas21[jxi]), imag(vec_fsas21[jxi])
	);
	fprintf(
	  p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	  vec_qschu1[jxi], vec_pschu1[jxi], vec_s0mag1[jxi]
	);
	fprintf(
	  p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
	  vec_cosav1[jxi], vec_raprs1[jxi]
	);
	fprintf(p_outfile, "  Fk=%15.7lE\n", vec_fk1[jxi]);
	if (inpol == 0) fprintf(p_outfile, "   LIN  1\n");
	else fprintf(p_outfile, "  CIRC  1\n");
	fprintf(p_outfile, " ----- SCS ----- ABS ----- EXS ----- ALBEDS --\n");
	fprintf(
	  p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
	  vec_scs2[jxi], vec_abs2[jxi], vec_exs2[jxi], vec_albeds2[jxi]
	);
	// fprintf(p_outfile, "INSERTION: SCASECM %5d%15.7E%15.7E%15.7E%15.7E\n",
	// 	1, alamb, vec_scs2[jxi], vec_abs2[jxi], vec_exs2[jxi]
	// );
	fprintf(p_outfile, " ---- SCS/GS -- ABC/GS -- EXS/GS ---\n");
	fprintf(
	  p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	  vec_scsrt2[jxi], vec_absrt2[jxi], vec_exsrt2[jxi]
	);
	fprintf(
	  p_outfile, "  FSAS(2,2)=%15.7lE%15.7lE   FSAS(1,2)=%15.7lE%15.7lE\n",
	  real(vec_fsas22[jxi]), imag(vec_fsas22[jxi]),
	  real(vec_fsas12[jxi]), imag(vec_fsas12[jxi])
	);
	fprintf(
	  p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	  vec_qschu2[jxi], vec_pschu2[jxi], vec_s0mag2[jxi]
	);
	fprintf(
	  p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
	  vec_cosav2[jxi], vec_raprs2[jxi]
	);
	fprintf(p_outfile, "  Fk=%15.7lE\n", vec_fk2[jxi]);
	fprintf(
	  p_outfile, "  (RE(FSAS(1,1))-RE(FSAS(2,2)))/RE(FSAS(1,1))=%15.7lE\n",
	  (real(vec_fsas11[jxi]) - real(vec_fsas22[jxi])) / real(vec_fsas11[jxi])
	);
	fprintf(
	  p_outfile, "  (IM(FSAS(1,1))-IM(FSAS(2,2)))/IM(FSAS(1,1))=%15.7lE\n",
	  (imag(vec_fsas11[jxi]) - imag(vec_fsas22[jxi])) / imag(vec_fsas11[jxi])
	);
	for (int jth = 0; jth < _num_theta; jth++) {
	  for (int jph = 0; jph < _num_phi; jph++) {
	    for (int jths = 0; jths < _num_thetas; jths++) {
	      for (int jphs = 0; jphs < _num_phis; jphs++) {
		bool goto290 = isam >= 0 && (jths > 0 || jphs > 0);
		int dir_index = ndirs * jxi + done_dirs;
		fprintf(
			p_outfile, "********** JTH =%3d, JPH =%3d, JTHS =%3d, JPHS =%3d ********************\n",
			jth + 1, jph + 1, jths + 1, jphs +1
		);
		fprintf(
		  p_outfile, "  TIDG=%10.3lE, PIDG=%10.3lE, TSDG=%10.3lE, PSDG=%10.3lE\n",
		  th + jth * thstp, ph + jph * phstp, ths + jths * thsstp, phs + jphs * phsstp
		);
		fprintf(p_outfile, "  SCAND=%10.3lE\n", vec_dir_scand[done_dirs]);
		fprintf(
		  p_outfile, "  CFMP=%15.7lE, SFMP=%15.7lE\n",
		  vec_dir_cfmp[done_dirs], vec_dir_sfmp[done_dirs]
		);
		fprintf(
		  p_outfile, "  CFSP=%15.7lE, SFSP=%15.7lE\n",
		  vec_dir_cfsp[done_dirs], vec_dir_sfsp[done_dirs]
		);
		if (isam >= 0) {
		  fprintf(
		    p_outfile, "  UNI=(%12.5lE,%12.5lE,%12.5lE)\n",
		    vec_dir_un[3 * done_dirs],
		    vec_dir_un[3 * done_dirs + 1],
		    vec_dir_un[3 * done_dirs + 2]
		  );
		  fprintf(
		    p_outfile, "  UNS=(%12.5lE,%12.5lE,%12.5lE)\n",
		    vec_dir_uns[3 * done_dirs],
		    vec_dir_uns[3 * done_dirs + 1],
		    vec_dir_uns[3 * done_dirs + 2]
		  );
		} else { // label 214
		  fprintf(
		    p_outfile, "  UN=(%12.5lE,%12.5lE,%12.5lE)\n\n",
		    vec_dir_un[3 * done_dirs],
		    vec_dir_un[3 * done_dirs + 1],
		    vec_dir_un[3 * done_dirs + 2]
		  );
		}
		fprintf(p_outfile, "     SINGLE SCATTERER\n");
		if (inpol == 0) {
		  fprintf(p_outfile, "   LIN -1\n");
		} else {
		  fprintf(p_outfile, "  CIRC -1\n");
		}
		// label 275
		fprintf(p_outfile, " ----- SCS ----- ABS ----- EXS ----- ALBEDS --\n");
		fprintf(
		  p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
		  vec_dir_scs1[dir_index],
		  vec_dir_abs1[dir_index],
		  vec_dir_exs1[dir_index],
		  vec_dir_albeds1[dir_index]
		);
		// fprintf(
		//   p_outfile, "INSERTION: SCASEC %5d%14.7lE%14.7lE%14.7lE%14.7lE\n",
		//   -1, alamb, vec_dir_scs1[dir_index], vec_dir_abs1[dir_index],
		//   vec_dir_exs1[dir_index]
		// );
		fprintf(p_outfile, " ---- SCS/GS -- ABS/GS -- EXS/GS ---\n");
		fprintf(
		  p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		  vec_dir_scsrt1[dir_index],
		  vec_dir_absrt1[dir_index],
		  vec_dir_exsrt1[dir_index]
		);
		fprintf(
		  p_outfile, "  FSAS(1,1)=%15.7lE%15.7lE   FSAS(2,1)=%15.7lE%15.7lE\n",
		  real(vec_dir_fsas11[dir_index]), imag(vec_dir_fsas11[dir_index]),
		  real(vec_dir_fsas21[dir_index]), imag(vec_dir_fsas21[dir_index])
		);
		fprintf(
		  p_outfile, "   SAS(1,1)=%15.7lE%15.7lE    SAS(2,1)=%15.7lE%15.7lE\n",
		  real(vec_dir_sas11[dir_index]), imag(vec_dir_sas11[dir_index]),
		  real(vec_dir_sas21[dir_index]), imag(vec_dir_sas21[dir_index])
		);
		fprintf(
		  p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
		  vec_dir_qschu1[dir_index],
		  vec_dir_pschu1[dir_index],
		  vec_dir_s0mag1[dir_index]
		);
		if (!goto290) {
		  fprintf(
		    p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
		    vec_dir_cosav1[dir_index], vec_dir_rapr1[dir_index]
		  );
		  fprintf(
		    p_outfile, "  Fl=%15.7lE, Fr=%15.7lE, Fk=%15.7lE\n",
		    vec_dir_fl1[dir_index], vec_dir_fr1[dir_index],
		    vec_dir_fk1[dir_index]
		  );
		  fprintf(
		    p_outfile, "  Fx=%15.7lE, Fy=%15.7lE, Fz=%15.7lE\n",
		    vec_dir_fx1[dir_index], vec_dir_fy1[dir_index],
		    vec_dir_fz1[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQEl=%15.7lE,  TQEr=%15.7lE,  TQEk=%15.7lE\n",
		    vec_dir_tqel1[dir_index], vec_dir_tqer1[dir_index],
		    vec_dir_tqek1[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQSl=%15.7lE,  TQSr=%15.7lE,  TQSk=%15.7lE\n",
		    vec_dir_tqsl1[dir_index], vec_dir_tqsr1[dir_index],
		    vec_dir_tqsk1[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQEx=%15.7lE,  TQEy=%15.7lE,  TQEz=%15.7lE\n",
		    vec_dir_tqex1[dir_index], vec_dir_tqey1[dir_index],
		    vec_dir_tqez1[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQSx=%15.7lE,  TQSy=%15.7lE,  TQSz=%15.7lE\n",
		    vec_dir_tqsx1[dir_index], vec_dir_tqsy1[dir_index],
		    vec_dir_tqsz1[dir_index]
		  );
		} // goto290 switch
		if (inpol == 0) {
		  fprintf(p_outfile, "   LIN  1\n");
		} else {
		  fprintf(p_outfile, "  CIRC  1\n");
		}
		// label 275
		fprintf(p_outfile, " ----- SCS ----- ABS ----- EXS ----- ALBEDS --\n");
		fprintf(
		  p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
		  vec_dir_scs2[dir_index],
		  vec_dir_abs2[dir_index],
		  vec_dir_exs2[dir_index],
		  vec_dir_albeds2[dir_index]
		);
		// fprintf(
		//   p_outfile, "INSERTION: SCASEC %5d%14.7lE%14.7lE%14.7lE%14.7lE\n",
		//   1, alamb, vec_dir_scs2[dir_index], vec_dir_abs2[dir_index],
		//   vec_dir_exs2[dir_index]
		// );
		fprintf(p_outfile, " ---- SCS/GS -- ABS/GS -- EXS/GS ---\n");
		fprintf(
		  p_outfile, " %14.7lE%15.7lE%15.7lE\n",
		  vec_dir_scsrt2[dir_index],
		  vec_dir_absrt2[dir_index],
		  vec_dir_exsrt2[dir_index]
		);
		fprintf(
		  p_outfile, "  FSAS(2,2)=%15.7lE%15.7lE   FSAS(1,2)=%15.7lE%15.7lE\n",
		  real(vec_dir_fsas22[dir_index]), imag(vec_dir_fsas22[dir_index]),
		  real(vec_dir_fsas12[dir_index]), imag(vec_dir_fsas12[dir_index])
		);
		fprintf(
		  p_outfile, "   SAS(2,2)=%15.7lE%15.7lE    SAS(1,2)=%15.7lE%15.7lE\n",
		  real(vec_dir_sas22[dir_index]), imag(vec_dir_sas22[dir_index]),
		  real(vec_dir_sas12[dir_index]), imag(vec_dir_sas12[dir_index])
		);
		fprintf(
		  p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
		  vec_dir_qschu2[dir_index],
		  vec_dir_pschu2[dir_index],
		  vec_dir_s0mag2[dir_index]
		);
		if (!goto290) {
		  fprintf(
		    p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
		    vec_dir_cosav2[dir_index], vec_dir_rapr2[dir_index]
		  );
		  fprintf(
		    p_outfile, "  Fl=%15.7lE, Fr=%15.7lE, Fk=%15.7lE\n",
		    vec_dir_fl2[dir_index], vec_dir_fr2[dir_index],
		    vec_dir_fk2[dir_index]
		  );
		  fprintf(
		    p_outfile, "  Fx=%15.7lE, Fy=%15.7lE, Fz=%15.7lE\n",
		    vec_dir_fx2[dir_index], vec_dir_fy2[dir_index],
		    vec_dir_fz2[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQEl=%15.7lE,  TQEr=%15.7lE,  TQEk=%15.7lE\n",
		    vec_dir_tqel2[dir_index], vec_dir_tqer2[dir_index],
		    vec_dir_tqek2[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQSl=%15.7lE,  TQSr=%15.7lE,  TQSk=%15.7lE\n",
		    vec_dir_tqsl2[dir_index], vec_dir_tqsr2[dir_index],
		    vec_dir_tqsk2[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQEx=%15.7lE,  TQEy=%15.7lE,  TQEz=%15.7lE\n",
		    vec_dir_tqex2[dir_index], vec_dir_tqey2[dir_index],
		    vec_dir_tqez2[dir_index]
		  );
		  fprintf(
		    p_outfile, "   TQSx=%15.7lE,  TQSy=%15.7lE,  TQSz=%15.7lE\n",
		    vec_dir_tqsx2[dir_index], vec_dir_tqsy2[dir_index],
		    vec_dir_tqsz2[dir_index]
		  );
		} // goto290 switch
		fprintf(
		  p_outfile, "  (RE(FSAS(1,1))-RE(FSAS(2,2)))/RE(FSAS(1,1))=%15.7lE\n",
		  (real(vec_dir_fsas11[dir_index]) - real(vec_dir_fsas22[dir_index])) / real(vec_dir_fsas11[dir_index])
		);
		fprintf(
		  p_outfile, "  (IM(FSAS(1,1))-IM(FSAS(2,2)))/IM(FSAS(1,1))=%15.7lE\n",
		  (imag(vec_dir_fsas11[dir_index]) - imag(vec_dir_fsas22[dir_index])) / imag(vec_dir_fsas11[dir_index])
		);
		fprintf(p_outfile, "  MULL\n");
		for (int i = 0; i < 4; i++) {
		  int mull_index = 16 * dir_index + 4 * i;
		  fprintf(
		    p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
		    vec_dir_mull[mull_index],
		    vec_dir_mull[mull_index + 1],
		    vec_dir_mull[mull_index + 2],
		    vec_dir_mull[mull_index + 3]
		  );
		} // i MULL loop
		fprintf(p_outfile, "  MULLLR\n");
		for (int i = 0; i < 4; i++) {
		  int mull_index = 16 * dir_index + 4 * i;
		  fprintf(
		    p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
		    vec_dir_mulllr[mull_index],
		    vec_dir_mulllr[mull_index + 1],
		    vec_dir_mulllr[mull_index + 2],
		    vec_dir_mulllr[mull_index + 3]
		  );
		} // i MULLR loop
		if (iavm != 0) {
		  fprintf(p_outfile, "     ENSEMBLE AVERAGE, MODE%2d\n", iavm);
		  if (inpol == 0) fprintf(p_outfile, "   LIN\n");
		  else fprintf(p_outfile, "  CIRC\n");
		  // label 318
		  fprintf(p_outfile, "  MULL\n");
		  for (int i = 0; i < 4; i++) {
		    int mul_dir_index = 16 * dir_index + 4 * i;
		    fprintf(
		      p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
		      vec_dir_mull[mul_dir_index],
		      vec_dir_mull[mul_dir_index + 1],
		      vec_dir_mull[mul_dir_index + 2],
		      vec_dir_mull[mul_dir_index + 3]
		    );
		  } // i MULL loop
		  fprintf(p_outfile, "  MULLLR\n");
		  for (int i = 0; i < 4; i++) {
		    int mul_dir_index = 16 * dir_index + 4 * i;
		    fprintf(
		      p_outfile, "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
		      vec_dir_mulllr[mul_dir_index],
		      vec_dir_mulllr[mul_dir_index + 1],
		      vec_dir_mulllr[mul_dir_index + 2],
		      vec_dir_mulllr[mul_dir_index + 3]
		    );
		  } // i MULLR loop
		} // iavm != 0 block
		done_dirs++;
	      } // jphs loop
	    } // jths loop
	  } // jph loop
	} // jth loop
    } // jxi wavelength loop
    fclose(p_outfile);
  } else {
    result = -1;
  }
  return result;
}

#ifdef MPI_VERSION
int InclusionOutputInfo::mpireceive(const mixMPI* mpidata, int pid) {
  int result = 0;
  int flag;
  int chk_nsph, chk_inpol, chk_iavm, chk_isam, chk_num_theta, chk_num_thetas;
  int chk_num_phi, chk_num_phis, chk_ndirs, chk_idfc, chk_configs;
  double chk_exri;
  MPI_Recv(&flag, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Proceed with the rest _only if_ flag==0, else nothing is to be received
  if (flag == 0) {
    MPI_Recv(&chk_nsph, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_inpol, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_iavm, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_isam, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_theta, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_thetas, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_phi, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_phis, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_ndirs, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_exri, 1, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_idfc, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_configs, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    result += (chk_nsph == nsph) ? 0 : 1;
    result += (chk_inpol == inpol) ? 0 : 1;
    result += (chk_iavm == iavm) ? 0 : 1;
    result += (chk_isam == isam) ? 0 : 1;
    result += (chk_num_theta == _num_theta) ? 0 : 1;
    result += (chk_num_thetas == _num_thetas) ? 0 : 1;
    result += (chk_num_phi == _num_phi) ? 0 : 1;
    result += (chk_num_phis == _num_phis) ? 0 : 1;
    result += (chk_ndirs == ndirs) ? 0 : 1;
    result += (chk_exri == exri) ? 0 : 1;
    result += (chk_idfc == idfc) ? 0 : 1;
    result += (chk_configs == configurations) ? 0 : 1;
    if (result == 0) {
      int xi1, offset, chunk_size;
      MPI_Send(&result, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD);
      MPI_Recv(&xi1, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Receive vectors of single values per scale
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = xi1 - _first_xi;
      MPI_Recv(vec_jxi + offset, chunk_size, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_ier + offset, chunk_size, MPI_SHORT, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_vk + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_xi + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_albeds1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_albeds2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scsrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scsrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_absrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_absrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exsrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exsrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschu1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschu2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschu1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschu2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0mag1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0mag2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_cosav1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_cosav2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_raprs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_raprs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fk1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fk2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsas11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsas21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsas22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsas12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors whose sizes depend on configurations
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * (configurations + 1);
      MPI_Recv(vec_sphere_sizes + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_ref_indices + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors whose sizes depend on directions and wavelengths
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * ndirs;
      MPI_Recv(vec_dir_scs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_scs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_abs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_abs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_exs1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_exs2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_albeds1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_albeds2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_scsrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_scsrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_absrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_absrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_exsrt1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_exsrt2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsas11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsas21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsas12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fsas22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qschu1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_qschu2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_pschu1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_pschu2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_s0mag1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_s0mag2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_cosav1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_cosav2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_rapr1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_rapr2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      MPI_Recv(vec_dir_fl1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fl2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fr1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fr2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fk1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fk2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fx1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fx2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fy1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fy2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fz1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fz2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqel1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqel2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqer1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqer2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqek1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqek2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqex1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqex2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqey1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqey2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqez1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqez2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsl1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsl2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsr1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsr2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsk1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsk2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsx1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsx2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsy1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsy2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsz1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_tqsz2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_mull + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_mulllr + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }  else {
      MPI_Send(&result, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD);
    }
  }
  return result;
}

int InclusionOutputInfo::mpisend(const mixMPI *mpidata) {
  int result = 0;
  int chunk_size;
  if (_skip_flag == 1) {
    // tell the receiver we are not sending anything
    int flag = 1;
    MPI_Send(&flag, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
  }
  else {
    // tell the receiver we are sending actual stuff
    int flag = 0;
    MPI_Send(&flag, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    // Send output metadata for configuration cross-check
    MPI_Send(&nsph, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&inpol, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&iavm, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&isam, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_theta, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_thetas, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_phi, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_phis, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&ndirs, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&exri, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&idfc, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&configurations, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    // Wait for process 0 to cross-check the configuration
    MPI_Recv(&result, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (result == 0) {
      // Process 0 confirmed the consistency of configuration. Send the data.
      // Send vectors of single values per scale
      MPI_Send(&_first_xi, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(&xi_block_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_jxi, xi_block_size, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_ier, xi_block_size, MPI_SHORT, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_vk, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_xi, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scs1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scs2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abs1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abs2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exs1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exs2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_albeds1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_albeds2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scsrt1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scsrt2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_absrt1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_absrt2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exsrt1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exsrt2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschu1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschu2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschu1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschu2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0mag1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0mag2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_cosav1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_cosav2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_raprs1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_raprs2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fk1, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fk2, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsas11, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsas21, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsas22, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsas12, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);

      // Send vectors whose sizes depend on configurations
      chunk_size = xi_block_size * (configurations + 1);
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_sizes, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_ref_indices, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);

      // Send vectors whose sizes depend on directions and wavelengths
      chunk_size = xi_block_size * ndirs;
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_scs1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_scs2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_abs1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_abs2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_exs1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_exs2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_albeds1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_albeds2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_scsrt1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_scsrt2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_absrt1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_absrt2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_exsrt1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_exsrt2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsas11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsas21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsas12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fsas22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qschu1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_qschu2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_pschu1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_pschu2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_s0mag1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_s0mag2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_cosav1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_cosav2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_rapr1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_rapr2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fl1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fl2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fr1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fr2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fk1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fk2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fx1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fx2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fy1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fy2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fz1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fz2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqel1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqel2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqer1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqer2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqek1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqek2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqex1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqex2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqey1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqey2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqez1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqez2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsl1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsl2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsr1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsr2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsk1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsk2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsx1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsx2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsy1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsy2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsz1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_tqsz2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_mull, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_mulllr, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    }
  }
  return result;
}
#endif // MPI_VERSION
// >>> END OF InclusionOutputInfo CLASS IMPLEMENTATION <<<

// >>> SphereOutputInfo CLASS IMPLEMENTATION <<<
SphereOutputInfo::SphereOutputInfo(
  ScattererConfiguration *sc, GeometryConfiguration *gc,
  const mixMPI *mpidata, int first_xi, int xi_length
) {
  _skip_flag = 0;
  _first_xi = first_xi;
  nsph = gc->number_of_spheres;
  lm = gc->l_max;
  inpol = gc->in_pol;
  npnt = gc->npnt;
  npntts = gc->npntts;
  isam = gc->isam;
  idfc = sc->idfc;
  th = gc->in_theta_start;
  thstp = gc->in_theta_step;
  thlst = gc->in_theta_end;
  _num_theta = (thstp == 0.0) ? 1 : 1 + (int)((thlst - th) / thstp);
  ths = gc->sc_theta_start;
  thsstp = gc->sc_theta_step;
  thslst = gc->sc_theta_end;
  _num_thetas = (thsstp == 0.0) ? 1 : 1 + (int)((thslst - ths) / thsstp);
  ph = gc->in_phi_start;
  phstp = gc->in_phi_step;
  phlst = gc->in_phi_end;
  _num_phi = (phstp == 0.0) ? 1 : 1 + (int)((phlst - ph) / phstp);
  phs = gc->sc_phi_start;
  phsstp = gc->sc_phi_step;
  phslst = gc->sc_phi_end;
  _num_phis = (phsstp == 0.0) ? 1 : 1 + (int)((phslst - phs) / phsstp);
  ndirs = _num_theta * _num_thetas * _num_phi * _num_phis;
  configurations = sc->configurations;
  double exdc = sc->exdc;
  exri = sqrt(exdc);
  nxi = sc->number_of_scales;
  xi_block_size = (xi_length == 0) ? nxi : xi_length;
  jwtm = gc->jwtm;
  lcalc = 0;
  arg = 0.0 + I * 0.0;
  vec_jxi = new int[xi_block_size]();
  vec_ier = new short[xi_block_size]();
  vec_vk = new double[xi_block_size]();
  vec_xi = new double[xi_block_size]();
  vec_sphere_sizes = new double[configurations * xi_block_size]();
  vec_sphere_ref_indices = new dcomplex[configurations * xi_block_size]();
  vec_scs = new double[configurations * xi_block_size]();
  vec_abs = new double[configurations * xi_block_size]();
  vec_exs = new double[configurations * xi_block_size]();
  vec_albeds = new double[configurations * xi_block_size]();
  vec_scsrt = new double[configurations * xi_block_size]();
  vec_absrt = new double[configurations * xi_block_size]();
  vec_exsrt = new double[configurations * xi_block_size]();
  vec_fsas = new dcomplex[configurations * xi_block_size]();
  vec_qschu = new double[configurations * xi_block_size]();
  vec_pschu = new double[configurations * xi_block_size]();
  vec_s0mag = new double[configurations * xi_block_size]();
  vec_cosav = new double[configurations * xi_block_size]();
  vec_raprs = new double[configurations * xi_block_size]();
  vec_tqek1 = new double[configurations * xi_block_size]();
  vec_tqek2 = new double[configurations * xi_block_size]();
  vec_tqsk1 = new double[configurations * xi_block_size]();
  vec_tqsk2 = new double[configurations * xi_block_size]();
  if (nsph == 1) {
    vec_fsat = NULL;
    vec_qschut = NULL;
    vec_pschut = NULL;
    vec_s0magt = NULL;
  } else {
    vec_fsat = new dcomplex[xi_block_size]();
    vec_qschut = new double[xi_block_size]();
    vec_pschut = new double[xi_block_size]();
    vec_s0magt = new double[xi_block_size]();
  }
  // Initialize directions (they are scale-independent)
  vec_dir_tidg = new double[_num_theta];
  vec_dir_pidg = new double[_num_phi];
  vec_dir_tsdg = new double[_num_thetas];
  vec_dir_psdg = new double[_num_phis];
  double cti = th, cpi = ph, cts = ths, cps = phs;
  for (int di = 0; di < _num_theta; di++) {
    vec_dir_tidg[di] = cti;
    cti += thstp;
  }
  for (int di = 0; di < _num_thetas; di++) {
    vec_dir_tsdg[di] = cts;
    cts += thsstp;
  }
  for (int di = 0; di < _num_phi; di++) {
    vec_dir_pidg[di] = cpi;
    cpi += phstp;
  }
  for (int di = 0; di < _num_phis; di++) {
    vec_dir_psdg[di] = cps;
    cps += phsstp;
  }
  vec_dir_scand = new double[ndirs]();
  vec_dir_cfmp = new double[ndirs]();
  vec_dir_cfsp = new double[ndirs]();
  vec_dir_sfmp = new double[ndirs]();
  vec_dir_sfsp = new double[ndirs]();
  vec_dir_un = new double[3 * ndirs]();
  vec_dir_uns = new double[3 * ndirs]();
  vec_dir_sas11 = new dcomplex[nsph * ndirs * xi_block_size]();
  vec_dir_sas21 = new dcomplex[nsph * ndirs * xi_block_size]();
  vec_dir_sas12 = new dcomplex[nsph * ndirs * xi_block_size]();
  vec_dir_sas22 = new dcomplex[nsph * ndirs * xi_block_size]();
  vec_dir_fx = new double[nsph * _num_theta * _num_phi * xi_block_size]();
  vec_dir_fy = new double[nsph * _num_theta * _num_phi * xi_block_size]();
  vec_dir_fz = new double[nsph * _num_theta * _num_phi * xi_block_size]();
  vec_dir_muls = new double[16 * nsph * ndirs * xi_block_size]();
  vec_dir_mulslr = new double[16 * nsph * ndirs * xi_block_size]();
}

SphereOutputInfo::SphereOutputInfo(const std::string &hdf5_name) {
  unsigned int flags = H5F_ACC_RDONLY;
  HDFFile *hdf_file = new HDFFile(hdf5_name, flags);
  herr_t status = hdf_file->get_status();
  string str_name, str_type;
  _skip_flag = 0;
  if (status == 0) {
    status = hdf_file->read("NSPH", "INT32_(1)", &nsph);
    status = hdf_file->read("LM", "INT32_(1)", &lm);
    status = hdf_file->read("INPOL", "INT32_(1)", &inpol);
    status = hdf_file->read("NPNT", "INT32_(1)", &npnt);
    status = hdf_file->read("NPNTTS", "INT32_(1)", &npntts);
    status = hdf_file->read("ISAM", "INT32_(1)", &isam);
    status = hdf_file->read("JWTM", "INT32_(1)", &jwtm);
    status = hdf_file->read("TH_START", "FLOAT64_(1)", &th);
    status = hdf_file->read("TH_STEP", "FLOAT64_(1)", &thstp);
    status = hdf_file->read("TH_END", "FLOAT64_(1)", &thlst);
    _num_theta = (thstp == 0.0) ? 1 : 1 + int((thlst - th) / thstp);
    status = hdf_file->read("THS_START", "FLOAT64_(1)", &ths);
    status = hdf_file->read("THS_STEP", "FLOAT64_(1)", &thsstp);
    status = hdf_file->read("THS_END", "FLOAT64_(1)", &thslst);
    _num_thetas = (thsstp == 0.0) ? 1 : 1 + int((thslst - ths) / thsstp);
    status = hdf_file->read("PH_START", "FLOAT64_(1)", &ph);
    status = hdf_file->read("PH_STEP", "FLOAT64_(1)", &phstp);
    status = hdf_file->read("PH_END", "FLOAT64_(1)", &phlst);
    _num_phi = (phstp == 0.0) ? 1 : 1 + int((phlst - ph) / phstp);
    status = hdf_file->read("PHS_START", "FLOAT64_(1)", &phs);
    status = hdf_file->read("PHS_STEP", "FLOAT64_(1)", &phsstp);
    status = hdf_file->read("PHS_END", "FLOAT64_(1)", &phslst);
    _num_phis = (phsstp == 0.0) ? 1 : 1 + int((phslst - phs) / phsstp);
    ndirs = _num_theta * _num_thetas * _num_phi * _num_phis;
    status = hdf_file->read("EXRI", "FLOAT64_(1)", &exri);
    status = hdf_file->read("NUM_CONF", "INT32_(1)", &configurations);
    status = hdf_file->read("IDFC", "INT32_(1)", &idfc);
    status = hdf_file->read("XI1", "INT32_(1)", &_first_xi);
    status = hdf_file->read("NXI", "INT32_(1)", &xi_block_size);
    nxi = (_first_xi == 1) ? xi_block_size : xi_block_size + _first_xi;
    lcalc = 0;
    arg = 0.0 + I * 0.0;
    str_type = "INT32_(" + to_string(xi_block_size) + ")";
    vec_jxi = new int[xi_block_size];
    status = hdf_file->read("VEC_JXI", str_type, vec_jxi);
    str_type = "INT16_(" + to_string(xi_block_size) + ")";
    vec_ier = new short[xi_block_size];
    status = hdf_file->read("VEC_IER", str_type, vec_ier);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    vec_vk = new double[xi_block_size];
    status = hdf_file->read("VEC_VK", str_type, vec_vk);
    vec_xi = new double[xi_block_size];
    status = hdf_file->read("VEC_XI", str_type, vec_xi);
    str_type = "FLOAT64_(" + to_string(configurations * xi_block_size) + ")";
    vec_sphere_sizes = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_SPH_SIZES", str_type, vec_sphere_sizes);
    str_type = "FLOAT64_(" + to_string(2 * configurations * xi_block_size) + ")";
    vec_sphere_ref_indices = new dcomplex[configurations * xi_block_size];
    status = hdf_file->read("VEC_SPH_REFRI", str_type, vec_sphere_ref_indices);
    str_type = "FLOAT64_(" + to_string(configurations * xi_block_size) + ")";
    vec_scs = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_SCS", str_type, vec_scs);
    vec_abs = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_ABS", str_type, vec_abs);
    vec_exs = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_EXS", str_type, vec_exs);
    vec_albeds = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_ALBEDS", str_type, vec_albeds);
    vec_scsrt = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_SCSRT", str_type, vec_scsrt);
    vec_absrt = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_ABSRT", str_type, vec_absrt);
    vec_exsrt = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_EXSRT", str_type, vec_exsrt);
    str_type = "FLOAT64_(" + to_string(2 * configurations * xi_block_size) + ")";
    vec_fsas = new dcomplex[configurations * xi_block_size];
    status = hdf_file->read("VEC_FSAS", str_type, vec_fsas);
    str_type = "FLOAT64_(" + to_string(configurations * xi_block_size) + ")";
    vec_qschu = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_QSCHU", str_type, vec_qschu);
    vec_pschu = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_PSCHU", str_type, vec_pschu);
    vec_s0mag = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_S0MAG", str_type, vec_s0mag);
    vec_cosav = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_COSAV", str_type, vec_cosav);
    vec_raprs = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_RAPRS", str_type, vec_raprs);
    vec_tqek1 = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_TQEK1", str_type, vec_tqek1);
    vec_tqek2 = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_TQEK2", str_type, vec_tqek2);
    vec_tqsk1 = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_TQSK1", str_type, vec_tqsk1);
    vec_tqsk2 = new double[configurations * xi_block_size];
    status = hdf_file->read("VEC_TQSK2", str_type, vec_tqsk2);
    if (nsph != 1) {
      str_type = "FLOAT64_(" + to_string(2 * xi_block_size) + ")";
      vec_fsat = new dcomplex[xi_block_size];
      status = hdf_file->read("VEC_FSAT", str_type, vec_fsat);
      str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
      vec_qschut = new double[xi_block_size];
      status = hdf_file->read("VEC_QSCHUT", str_type, vec_qschut);
      vec_pschut = new double[xi_block_size];
      status = hdf_file->read("VEC_PSCHUT", str_type, vec_pschut);
      vec_s0magt = new double[xi_block_size];
      status = hdf_file->read("VEC_S0MAGT", str_type, vec_s0magt);
    } else {
      vec_fsat = NULL;
      vec_qschut = NULL;
      vec_pschut = NULL;
      vec_s0magt = NULL;
    }
    // Initialize directions (they are scale-independent)
    vec_dir_tidg = new double[_num_theta];
    vec_dir_tsdg = new double[_num_thetas];
    vec_dir_pidg = new double[_num_phi];
    vec_dir_psdg = new double[_num_phis];
    double cti = th, cpi = ph, cts = ths, cps = phs;
    for (int di = 0; di < _num_theta; di++) {
      vec_dir_tidg[di] = cti;
      cti += thstp;
    }
    for (int di = 0; di < _num_thetas; di++) {
      vec_dir_tsdg[di] = cts;
      cts += thsstp;
    }
    for (int di = 0; di < _num_phi; di++) {
      vec_dir_pidg[di] = cpi;
      cpi += phstp;
    }
    for (int di = 0; di < _num_phis; di++) {
      vec_dir_psdg[di] = cps;
      cps += phsstp;
    }
    str_type = "FLOAT64_(" + to_string(ndirs) + ")";
    vec_dir_scand = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SCAND", str_type, vec_dir_scand);
    vec_dir_cfmp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_CFMP", str_type, vec_dir_cfmp);
    vec_dir_sfmp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SFMP", str_type, vec_dir_sfmp);
    vec_dir_cfsp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_CFSP", str_type, vec_dir_cfsp);
    vec_dir_sfsp = new double[ndirs];
    status = hdf_file->read("VEC_DIR_SFSP", str_type, vec_dir_sfsp);
    str_type = "FLOAT64_(" + to_string(3 * ndirs) + ")";
    vec_dir_un = new double[3 * ndirs];
    status = hdf_file->read("VEC_DIR_UN", str_type, vec_dir_un);
    vec_dir_uns = new double[3 * ndirs];
    status = hdf_file->read("VEC_DIR_UNS", str_type, vec_dir_uns);
    str_type = "FLOAT64_(" + to_string(2 * nsph * ndirs * xi_block_size) + ")";
    vec_dir_sas11 = new dcomplex[nsph * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS11", str_type, vec_dir_sas11);
    vec_dir_sas21 = new dcomplex[nsph * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS21", str_type, vec_dir_sas21);
    vec_dir_sas12 = new dcomplex[nsph * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS12", str_type, vec_dir_sas12);
    vec_dir_sas22 = new dcomplex[nsph * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_SAS22", str_type, vec_dir_sas22);
    str_type = "FLOAT64_(" + to_string(nsph * _num_theta * _num_phi * xi_block_size) + ")";
    vec_dir_fx = new double[nsph * _num_theta * _num_phi * xi_block_size];
    status = hdf_file->read("VEC_DIR_FX", str_type, vec_dir_fx);
    vec_dir_fy = new double[nsph * _num_theta * _num_phi * xi_block_size];
    status = hdf_file->read("VEC_DIR_FY", str_type, vec_dir_fy);
    vec_dir_fz = new double[nsph * _num_theta * _num_phi * xi_block_size];
    status = hdf_file->read("VEC_DIR_FZ", str_type, vec_dir_fz);
    str_type = "FLOAT64_(" + to_string(16 * nsph * ndirs * xi_block_size) + ")";
    vec_dir_muls = new double[16 * nsph * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULS", str_type, vec_dir_muls);
    vec_dir_mulslr = new double[16 * nsph * ndirs * xi_block_size];
    status = hdf_file->read("VEC_DIR_MULSLR", str_type, vec_dir_mulslr);
    status = hdf_file->close();
    delete hdf_file;
  } else {
    if (hdf_file != NULL) delete hdf_file;
    UnrecognizedFormatException ex("Error: " + hdf5_name + " not recognized as a valid HDF5 file!");
    throw ex;
  }
}

SphereOutputInfo::SphereOutputInfo(const int flag) {
  /*
    create a dummy placeholder just to know I should skip MPI_Send and MPI_Recv
  */
  if (flag == 1) {
    _skip_flag = 1;
  } else {
    UnrecognizedOutputInfo ex(flag);
    throw ex;
  }
}

SphereOutputInfo::~SphereOutputInfo() {
  if (_skip_flag != 1) {
    delete[] vec_jxi;
    delete[] vec_ier;
    delete[] vec_vk;
    delete[] vec_xi;
    delete[] vec_sphere_sizes;
    delete[] vec_sphere_ref_indices;
    delete[] vec_scs;
    delete[] vec_abs;
    delete[] vec_exs;
    delete[] vec_albeds;
    delete[] vec_scsrt;
    delete[] vec_absrt;
    delete[] vec_exsrt;
    delete[] vec_fsas;
    delete[] vec_qschu;
    delete[] vec_pschu;
    delete[] vec_s0mag;
    delete[] vec_cosav;
    delete[] vec_raprs;
    delete[] vec_tqek1;
    delete[] vec_tqek2;
    delete[] vec_tqsk1;
    delete[] vec_tqsk2;
    if (nsph != 1) {
      delete[] vec_fsat;
      delete[] vec_qschut;
      delete[] vec_pschut;
      delete[] vec_s0magt;
    }
    delete[] vec_dir_tidg;
    delete[] vec_dir_pidg;
    delete[] vec_dir_tsdg;
    delete[] vec_dir_psdg;
    delete[] vec_dir_scand;
    delete[] vec_dir_cfmp;
    delete[] vec_dir_cfsp;
    delete[] vec_dir_sfmp;
    delete[] vec_dir_sfsp;
    delete[] vec_dir_un;
    delete[] vec_dir_uns;
    delete[] vec_dir_sas11;
    delete[] vec_dir_sas21;
    delete[] vec_dir_sas12;
    delete[] vec_dir_sas22;
    delete[] vec_dir_fx;
    delete[] vec_dir_fy;
    delete[] vec_dir_fz;
    delete[] vec_dir_muls;
    delete[] vec_dir_mulslr;
  }
}

long SphereOutputInfo::compute_size(
  ScattererConfiguration *sc, GeometryConfiguration *gc,
  int first_xi, int xi_length
) {
  // Size of the configuration set
  long result = 18 * sizeof(int);
  result += 12 * sizeof(double);
  result += 47 * sizeof(long);
  result += sizeof(dcomplex);
  // Get configuration parameters
  int _nsph = gc->number_of_spheres;
  double _th = gc->in_theta_start;
  double _thstp = gc->in_theta_step;
  double _thlst = gc->in_theta_end;
  int num_theta = (_thstp == 0.0) ? 1 : 1 + (int)((_thlst - _th) / _thstp);
  double _ths = gc->sc_theta_start;
  double _thsstp = gc->sc_theta_step;
  double _thslst = gc->sc_theta_end;
  int num_thetas = (_thsstp == 0.0) ? 1 : 1 + (int)((_thslst - _ths) / _thsstp);
  double _ph = gc->in_phi_start;
  double _phstp = gc->in_phi_step;
  double _phlst = gc->in_phi_end;
  int num_phi = (_phstp == 0.0) ? 1 : 1 + (int)((_phlst - _ph) / _phstp);
  double _phs = gc->sc_phi_start;
  double _phsstp = gc->sc_phi_step;
  double _phslst = gc->sc_phi_end;
  int num_phis = (_phsstp == 0.0) ? 1 : 1 + (int)((_phslst - _phs) / _phsstp);
  int _ndirs = num_theta * num_thetas * num_phi * num_phis;
  int _nxi = sc->number_of_scales;
  int _xi_block_size = (xi_length == 0) ? _nxi : xi_length;
  int _nconf = sc->configurations;
  // Size of the data set
  result += _xi_block_size * (sizeof(short) + sizeof(int));
  result += 2 * _xi_block_size * sizeof(double);
  result += 16 * _nconf * _xi_block_size * sizeof(double);
  result += 2 * _nconf * _xi_block_size * sizeof(dcomplex);
  result += (num_theta + num_thetas + num_phi + num_phis) * sizeof(double);
  result += 11 * _ndirs * sizeof(double);
  result += 4 * _nsph * _ndirs * _xi_block_size * sizeof(dcomplex);
  result += 3 * _nsph * num_theta * num_phi * sizeof(double);
  result += 32 * _nsph * _ndirs * _xi_block_size * sizeof(double);
  if (_nsph != 1) {
    result += _xi_block_size * sizeof(dcomplex);
    result += 3 * _xi_block_size * sizeof(double);
  }
  return result;
}

long SphereOutputInfo::compute_size() {
  // Size of the configuration set
  long result = 18 * sizeof(int);
  result += 12 * sizeof(double);
  result += 47 * sizeof(long);
  result += sizeof(dcomplex);
  // Size of the data set
  result += xi_block_size * (sizeof(short) + sizeof(int));
  result += 2 * xi_block_size * sizeof(double);
  result += 16 * configurations * xi_block_size * sizeof(double);
  result += 2 * configurations * xi_block_size * sizeof(dcomplex);
  result += (_num_theta + _num_thetas + _num_phi + _num_phis) * sizeof(double);
  result += 11 * ndirs * sizeof(double);
  result += 4 * nsph * ndirs * xi_block_size * sizeof(dcomplex);
  result += 3 * nsph * _num_theta * _num_phi * sizeof(double);
  result += 32 * nsph * ndirs * xi_block_size * sizeof(double);
  if (nsph != 1) {
    result += xi_block_size * sizeof(dcomplex);
    result += 3 * xi_block_size * sizeof(double);
  }
  return result;
}

int SphereOutputInfo::insert(const SphereOutputInfo &rhs) {
  int result = 0;
  if (rhs.skip_flag != 1) {
    result += (rhs.nsph == nsph) ? 0 : 1;
    result += (rhs.inpol == inpol) ? 0 : 1;
    result += (rhs.isam == isam) ? 0 : 1;
    result += (rhs._num_theta == _num_theta) ? 0 : 1;
    result += (rhs._num_thetas == _num_thetas) ? 0 : 1;
    result += (rhs._num_phi == _num_phi) ? 0 : 1;
    result += (rhs._num_phis == _num_phis) ? 0 : 1;
    result += (rhs.ndirs == ndirs) ? 0 : 1;
    result += (rhs.exri == exri) ? 0 : 1;
    result += (rhs.idfc == idfc) ? 0 : 1;
    result += (rhs.configurations == configurations) ? 0 : 1;
    if (result == 0) {
      int offset, chunk_size, xi1;
      xi1 = rhs._first_xi;
      // Insert vectors whose sizes depend on wavelengths
      offset = xi1 - _first_xi;
      chunk_size = rhs.xi_block_size;
      memcpy(vec_jxi + offset, rhs.vec_jxi, chunk_size * sizeof(int));
      memcpy(vec_ier + offset, rhs.vec_ier, chunk_size * sizeof(short));
      memcpy(vec_vk + offset, rhs.vec_vk, chunk_size * sizeof(double));
      memcpy(vec_xi + offset, rhs.vec_xi, chunk_size * sizeof(double));
      if (nsph != 1) {
	memcpy(vec_fsat + offset, rhs.vec_fsat, chunk_size * sizeof(dcomplex));
	memcpy(vec_qschut + offset, rhs.vec_qschut, chunk_size * sizeof(double));
	memcpy(vec_pschut + offset, rhs.vec_pschut, chunk_size * sizeof(double));
	memcpy(vec_s0magt + offset, rhs.vec_s0magt, chunk_size * sizeof(double));
      }

      // Insert vectors whose sizes depend on configurations and wavelengths
      offset = (xi1 - _first_xi) * configurations;
      chunk_size = rhs.xi_block_size * configurations;
      memcpy(vec_sphere_sizes + offset, rhs.vec_sphere_sizes, chunk_size * sizeof(double));
      memcpy(vec_sphere_ref_indices + offset, rhs.vec_sphere_ref_indices, chunk_size * sizeof(dcomplex));
      memcpy(vec_scs + offset, rhs.vec_scs, chunk_size * sizeof(double));
      memcpy(vec_abs + offset, rhs.vec_abs, chunk_size * sizeof(double));
      memcpy(vec_exs + offset, rhs.vec_exs, chunk_size * sizeof(double));
      memcpy(vec_albeds + offset, rhs.vec_albeds, chunk_size * sizeof(double));
      memcpy(vec_scsrt + offset, rhs.vec_scsrt, chunk_size * sizeof(double));
      memcpy(vec_absrt + offset, rhs.vec_absrt, chunk_size * sizeof(double));
      memcpy(vec_exsrt + offset, rhs.vec_exsrt, chunk_size * sizeof(double));
      memcpy(vec_fsas + offset, rhs.vec_fsas, chunk_size * sizeof(dcomplex));
      memcpy(vec_qschu + offset, rhs.vec_qschu, chunk_size * sizeof(double));
      memcpy(vec_pschu + offset, rhs.vec_pschu, chunk_size * sizeof(double));
      memcpy(vec_s0mag + offset, rhs.vec_s0mag, chunk_size * sizeof(double));
      memcpy(vec_cosav + offset, rhs.vec_cosav, chunk_size * sizeof(double));
      memcpy(vec_raprs + offset, rhs.vec_raprs, chunk_size * sizeof(double));
      memcpy(vec_tqek1 + offset, rhs.vec_tqek1, chunk_size * sizeof(double));
      memcpy(vec_tqek2 + offset, rhs.vec_tqek2, chunk_size * sizeof(double));
      memcpy(vec_tqsk1 + offset, rhs.vec_tqsk1, chunk_size * sizeof(double));
      memcpy(vec_tqsk2 + offset, rhs.vec_tqsk2, chunk_size * sizeof(double));
    
      // Insert vectors whose sizes depend on NSPH, directions and wavelengths
      offset = (xi1 - _first_xi) * nsph * ndirs;
      chunk_size = rhs.xi_block_size * nsph * ndirs;
      memcpy(vec_dir_sas11 + offset, rhs.vec_dir_sas11, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas21 + offset, rhs.vec_dir_sas21, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas12 + offset, rhs.vec_dir_sas12, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_sas22 + offset, rhs.vec_dir_sas22, chunk_size * sizeof(dcomplex));
      memcpy(vec_dir_muls + 16 * offset, rhs.vec_dir_muls, 16 * chunk_size * sizeof(double));
      memcpy(vec_dir_mulslr + 16 * offset, rhs.vec_dir_mulslr, 16 * chunk_size * sizeof(double));

      // Insert vectors whose sizes depend on NSPH, incidence directions and wavelengths
      offset = (xi1 - _first_xi) * nsph * _num_theta * _num_phi;
      chunk_size = rhs.xi_block_size  * nsph * _num_theta * _num_phi;
      memcpy(vec_dir_fx + offset, rhs.vec_dir_fx, chunk_size * sizeof(double));
      memcpy(vec_dir_fy + offset, rhs.vec_dir_fy, chunk_size * sizeof(double));
      memcpy(vec_dir_fz + offset, rhs.vec_dir_fz, chunk_size * sizeof(double));
      // TODO: fix the vector sizes in HDF5 writer and MPI communicators
    }
  }
  return result;
}

int SphereOutputInfo::write(const std::string &output, const std::string &format) {
  int result = 0;
  if (format.compare("LEGACY") == 0) {
    result = write_legacy(output);
  } else if (format.compare("HDF5") == 0) {
    result = write_hdf5(output);
  } else {
    string message = "Unknown format mode: \"" + format + "\"";
    throw UnrecognizedConfigurationException(message);
  }
  return result;
}

int SphereOutputInfo::write_hdf5(const std::string &file_name) {
  List<string> *rec_name_list = new List<string>(1);
  List<string> *rec_type_list = new List<string>(1);
  List<void *> *rec_ptr_list = new List<void *>(1);
  string str_type, str_name;
  rec_name_list->set(0, "NSPH");
  rec_type_list->set(0, "INT32_(1)");
  rec_ptr_list->set(0, &nsph);
  rec_name_list->append("LM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&lm);
  rec_name_list->append("INPOL");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&inpol);
  rec_name_list->append("NPNT");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&npnt);
  rec_name_list->append("NPNTTS");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&npntts);
  rec_name_list->append("ISAM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&isam);
  rec_name_list->append("JWTM");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&jwtm);
  rec_name_list->append("TH_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&th);
  rec_name_list->append("TH_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thstp);
  rec_name_list->append("TH_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thlst);
  rec_name_list->append("THS_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&ths);
  rec_name_list->append("THS_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thsstp);
  rec_name_list->append("THS_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&thslst);
  rec_name_list->append("PH_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&ph);
  rec_name_list->append("PH_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phstp);
  rec_name_list->append("PH_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phlst);
  rec_name_list->append("PHS_START");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phs);
  rec_name_list->append("PHS_STEP");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phsstp);
  rec_name_list->append("PHS_END");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&phslst);
  rec_name_list->append("EXRI");
  rec_type_list->append("FLOAT64_(1)");
  rec_ptr_list->append(&exri);
  rec_name_list->append("IDFC");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&idfc);
  rec_name_list->append("NUM_CONF");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&configurations);
  rec_name_list->append("XI1");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&_first_xi);
  rec_name_list->append("NXI");
  rec_type_list->append("INT32_(1)");
  rec_ptr_list->append(&xi_block_size);
  rec_name_list->append("VEC_JXI");
  rec_type_list->append("INT32_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_jxi);
  rec_name_list->append("VEC_IER");
  rec_type_list->append("INT16_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_ier);
  rec_name_list->append("VEC_VK");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_vk);
  rec_name_list->append("VEC_XI");
  rec_type_list->append("FLOAT64_(" + to_string(xi_block_size) + ")");
  rec_ptr_list->append(vec_xi);
  rec_name_list->append("VEC_SPH_SIZES");
  rec_type_list->append("FLOAT64_(" + to_string(configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_sphere_sizes);
  rec_name_list->append("VEC_SPH_REFRI");
  rec_type_list->append("FLOAT64_(" + to_string(2 * configurations * xi_block_size) + ")");
  rec_ptr_list->append(vec_sphere_ref_indices);
  str_type = "FLOAT64_(" + to_string(configurations * xi_block_size) + ")";
  rec_name_list->append("VEC_SCS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_scs);
  rec_name_list->append("VEC_ABS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_abs);
  rec_name_list->append("VEC_EXS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_exs);
  rec_name_list->append("VEC_ALBEDS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_albeds);
  rec_name_list->append("VEC_SCSRT");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_scsrt);
  rec_name_list->append("VEC_ABSRT");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_absrt);
  rec_name_list->append("VEC_EXSRT");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_exsrt);
  str_type = "FLOAT64_(" + to_string(2 * configurations * xi_block_size) + ")";
  rec_name_list->append("VEC_FSAS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_fsas);
  str_type = "FLOAT64_(" + to_string(configurations * xi_block_size) + ")";
  rec_name_list->append("VEC_QSCHU");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_qschu);
  rec_name_list->append("VEC_PSCHU");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_pschu);
  rec_name_list->append("VEC_S0MAG");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_s0mag);
  rec_name_list->append("VEC_COSAV");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_cosav);
  rec_name_list->append("VEC_RAPRS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_raprs);
  rec_name_list->append("VEC_TQEK1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_tqek1);
  rec_name_list->append("VEC_TQEK2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_tqek2);
  rec_name_list->append("VEC_TQSK1");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_tqsk1);
  rec_name_list->append("VEC_TQSK2");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_tqsk2);
  if (nsph != 1) {
    str_type = "FLOAT64_(" + to_string(2 * xi_block_size) + ")";
    rec_name_list->append("VEC_FSAT");
    rec_type_list->append(str_type);
    rec_ptr_list->append(vec_fsat);
    str_type = "FLOAT64_(" + to_string(xi_block_size) + ")";
    rec_name_list->append("VEC_QSCHUT");
    rec_type_list->append(str_type);
    rec_ptr_list->append(vec_qschut);
    rec_name_list->append("VEC_PSCHUT");
    rec_type_list->append(str_type);
    rec_ptr_list->append(vec_pschut);
    rec_name_list->append("VEC_S0MAGT");
    rec_type_list->append(str_type);
    rec_ptr_list->append(vec_s0magt);
  }
  str_type = "FLOAT64_(" + to_string(ndirs) + ")";
  rec_name_list->append("VEC_DIR_SCAND");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_scand);
  rec_name_list->append("VEC_DIR_CFMP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_cfmp);
  rec_name_list->append("VEC_DIR_CFSP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_cfsp);
  rec_name_list->append("VEC_DIR_SFMP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sfmp);
  rec_name_list->append("VEC_DIR_SFSP");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sfsp);
  str_type = "FLOAT64_(" + to_string(3 * ndirs) + ")";
  rec_name_list->append("VEC_DIR_UN");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_un);
  rec_name_list->append("VEC_DIR_UNS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_uns);
  str_type = "FLOAT64_(" + to_string(2 * nsph * ndirs * xi_block_size) + ")";
  rec_name_list->append("VEC_DIR_SAS11");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas11);
  rec_name_list->append("VEC_DIR_SAS21");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas21);
  rec_name_list->append("VEC_DIR_SAS12");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas12);
  rec_name_list->append("VEC_DIR_SAS22");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_sas22);
  str_type = "FLOAT64_(" + to_string(nsph * _num_theta * _num_phi * xi_block_size) + ")";
  rec_name_list->append("VEC_DIR_FX");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fx);
  rec_name_list->append("VEC_DIR_FY");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fy);
  rec_name_list->append("VEC_DIR_FZ");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_fz);
  str_type = "FLOAT64_(" + to_string(16 * nsph * ndirs * xi_block_size) + ")";
  rec_name_list->append("VEC_DIR_MULS");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_muls);
  rec_name_list->append("VEC_DIR_MULSLR");
  rec_type_list->append(str_type);
  rec_ptr_list->append(vec_dir_mulslr);
  
  // Convert the lists to arrays and write them to HDF5
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
  delete[] rec_names;
  delete[] rec_types;
  delete[] rec_pointers;
  delete schema;
  delete hdf_file;
  return 0;
}

int SphereOutputInfo::write_legacy(const std::string &file_name) {
  const dcomplex cc0 = 0.0 + I * 0.0;
  int result = 0;
  int nks = _num_thetas * _num_phis;
  FILE *p_outfile = fopen(file_name.c_str(), "w");
  if (p_outfile != NULL) {
    if (vec_jxi[0] == 1) {
      fprintf(p_outfile, " READ(IR,*)NSPH,LM,INPOL,NPNT,NPNTTS,ISAM\n");
      fprintf(
        p_outfile, " %5d%5d%5d%5d%5d%5d\n",
        nsph, lm, inpol, npnt, npntts, isam
      );
      fprintf(p_outfile, " READ(IR,*)TH,THSTP,THLST,THS,THSSTP,THSLST\n");
      fprintf(
        p_outfile, "  %9.3lE %9.3lE %9.3lE %9.3lE %9.3lE %9.3lE\n",
        th, thstp, thlst, ths, thsstp, thslst
      );
      fprintf(p_outfile, " READ(IR,*)PH,PHSTP,PHLST,PHS,PHSSTP,PHSLST\n");
      fprintf(
        p_outfile, "  %9.3lE %9.3lE %9.3lE %9.3lE %9.3lE %9.3lE\n",
        ph, phstp, phlst, phs, phsstp, phslst
      );
      fprintf(p_outfile, " READ(IR,*)JWTM\n");
      fprintf(p_outfile, " %5d\n", jwtm);
      fprintf(p_outfile, "  READ(ITIN)NSPHT\n");
      fprintf(p_outfile, "  READ(ITIN)(IOG(I),I=1,NSPH)\n");
      fprintf(p_outfile, "  READ(ITIN)EXDC,WP,XIP,IDFC,NXI\n");
      fprintf(p_outfile, "  READ(ITIN)(XIV(I),I=1,NXI)\n");
      fprintf(p_outfile, "  READ(ITIN)NSHL(I),ROS(I)\n");
      fprintf(p_outfile, "  READ(ITIN)(RCF(I,NS),NS=1,NSH)\n \n");
      fprintf(p_outfile, "  REFR. INDEX OF EXTERNAL MEDIUM=%15.7lE\n", exri);
      if (inpol == 0) fprintf(p_outfile, "   LIN\n \n");
      else fprintf(p_outfile, "  CIRC\n \n");
      if (idfc < 0) {
	fprintf(p_outfile, "  VK=%15.7lE, XI IS SCALE FACTOR FOR LENGTHS\n \n", vec_vk[0]);
      }
    } // End preamble writing
    // Wavelength loop
    for (int jxi = 0; jxi < xi_block_size; jxi++) {
      int done_dirs = 0;
      fprintf(p_outfile, "========== JXI =%3d ====================\n", jxi + 1);
      if (idfc >= 0) {
	fprintf(p_outfile, "  VK=%15.7lE, XI=%15.7lE\n", vec_vk[jxi], vec_xi[jxi]);
      } else { // IDFC < 0
	fprintf(p_outfile, "  XI=%15.7lE\n", vec_xi[jxi]);
      }
      if (vec_ier[jxi] == 1) {
	fprintf(p_outfile, "  STOP IN DME\n");
	fprintf(
	  p_outfile, "  AT %1d LCALC=%3d TOO SMALL WITH ARG=%15.7lE+i(%15.7lE)\n",
	  (int)vec_ier[jxi], lcalc, real(arg), imag(arg)
	);
      }
      for (int ci = 0; ci < configurations; ci++) {
	int cindex = jxi * configurations + ci;
	fprintf(p_outfile, "     SPHERE %2d\n", ci + 1);
	if (vec_sphere_ref_indices[cindex] == cc0) {
	  fprintf(p_outfile, "  SIZE=%15.7lE\n", vec_sphere_sizes[cindex]);
	} else {
	  fprintf(
		  p_outfile, "  SIZE=%15.7lE, REFRACTIVE INDEX=%15.7lE%15.7lE\n",
		  vec_sphere_sizes[cindex], real(vec_sphere_ref_indices[cindex]),
		  imag(vec_sphere_ref_indices[cindex])
	  );
	}
	fprintf(p_outfile, " ----- SCS ----- ABS ----- EXS ----- ALBEDS --\n");
	fprintf(
	  p_outfile, " %14.7lE%15.7lE%15.7lE%15.7lE\n",
	  vec_scs[cindex], vec_abs[cindex], vec_exs[cindex], vec_albeds[cindex]
        );
	fprintf(p_outfile, " ---- SCS/GS -- ABS/GS -- EXS/GS ---\n");
	fprintf(
	  p_outfile, " %14.7lE%15.7lE%15.7lE\n",
	  vec_scsrt[cindex], vec_absrt[cindex], vec_exsrt[cindex]
        );
	fprintf(
          p_outfile, "  FSAS=%15.7lE%15.7lE\n",
	  real(vec_fsas[cindex]), imag(vec_fsas[cindex])
        );
	fprintf(
	  p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	  vec_qschu[cindex], vec_pschu[cindex], vec_s0mag[cindex]
        );
	fprintf(
	  p_outfile, "  COSAV=%15.7lE, RAPRS=%15.7lE\n",
	  vec_cosav[cindex], vec_raprs[cindex]
        );
	fprintf(
	  p_outfile, "  IPO= 1, TQEk=%15.7lE, TQSk=%15.7lE\n",
	  vec_tqek1[cindex], vec_tqsk1[cindex]
	);
	fprintf(
	  p_outfile, "  IPO= 2, TQEk=%15.7lE, TQSk=%15.7lE\n",
	  vec_tqek2[cindex], vec_tqsk2[cindex]
        );
      } // ci configuration loop
      if (nsph != 1) {
	fprintf(
	  p_outfile, "  FSAT=(%15.7lE,%15.7lE)\n",
	  real(vec_fsat[jxi]), imag(vec_fsat[jxi])
	);
	fprintf(
	  p_outfile, "  QSCHU=%15.7lE, PSCHU=%15.7lE, S0MAG=%15.7lE\n",
	  vec_qschut[jxi], vec_pschut[jxi], vec_s0magt[jxi]
	);
      }
      for (int jth = 0; jth < _num_theta; jth++) {
	for (int jph = 0; jph < _num_phi; jph++) {
	  for (int jths = 0; jths < _num_thetas; jths++) {
	    for (int jphs = 0; jphs < _num_phis; jphs++) {
	      int dir_index = ndirs * jxi + done_dirs;
	      bool goto190 = (nks == 1) && ((jxi > 0) || (jth > 0) || (jph > 0));
	      fprintf(
		p_outfile, "********** JTH =%3d, JPH =%3d, JTHS =%3d, JPHS =%3d ********************\n",
		jth + 1, jph + 1, jths + 1, jphs + 1
	      );
	      fprintf(
		p_outfile, "  TIDG=%10.3lE, PIDG=%10.3lE, TSDG=%10.3lE, PSDG=%10.3lE\n",
		th + jth * thstp,
		ph + jph * phstp,
		ths + jths * thsstp,
		phs + jphs * phsstp
	      );
	      fprintf(p_outfile, "  SCAND=%10.3lE\n", vec_dir_scand[done_dirs]);
	      fprintf(
		p_outfile, "  CFMP=%15.7lE, SFMP=%15.7lE\n",
		vec_dir_cfmp[done_dirs], vec_dir_sfmp[done_dirs]
	      );
	      fprintf(
		p_outfile, "  CFSP=%15.7lE, SFSP=%15.7lE\n",
		vec_dir_cfsp[done_dirs], vec_dir_sfsp[done_dirs]
	      );
	      if (isam >= 0) {
		fprintf(
		  p_outfile, "  UNI=(%12.5lE,%12.5lE,%12.5lE)\n",
		  vec_dir_un[3 * done_dirs],
		  vec_dir_un[3 * done_dirs + 1],
		  vec_dir_un[3 * done_dirs + 2]
		);
		fprintf(
		  p_outfile, "  UNS=(%12.5lE,%12.5lE,%12.5lE)\n",
		  vec_dir_uns[3 * done_dirs],
		  vec_dir_uns[3 * done_dirs + 1],
		  vec_dir_uns[3 * done_dirs + 2]
		);
	      } else {
		fprintf(
		  p_outfile, "  UN=(%12.5lE,%12.5lE,%12.5lE)\n",
		  vec_dir_un[3 * done_dirs],
		  vec_dir_un[3 * done_dirs + 1],
		  vec_dir_un[3 * done_dirs + 2]
		);
	      }
	      for (int i = 0; i < nsph; i++) {
		int cindex = jxi * nsph * ndirs + nsph * done_dirs + i;
		fprintf(p_outfile, "     SPHERE %2d\n", i + 1);
		fprintf(
		  p_outfile, "  SAS(1,1)=%15.7lE%15.7lE, SAS(2,1)=%15.7lE%15.7lE\n",
		  real(vec_dir_sas11[cindex]),
		  imag(vec_dir_sas11[cindex]),
		  real(vec_dir_sas21[cindex]),
		  imag(vec_dir_sas21[cindex])
		);
		fprintf(
		  p_outfile, "  SAS(1,2)=%15.7lE%15.7lE, SAS(2,2)=%15.7lE%15.7lE\n",
		  real(vec_dir_sas12[cindex]),
		  imag(vec_dir_sas12[cindex]),
		  real(vec_dir_sas22[cindex]),
		  imag(vec_dir_sas22[cindex])
		);
		if (jths == 0 && jphs == 0) {
		  fprintf(
		    p_outfile, "  Fx=%15.7lE, Fy=%15.7lE, Fz=%15.7lE\n",
		    vec_dir_fx[jxi * nsph * _num_theta * _num_phi + jth * nsph * _num_phi + jph * nsph + i],
		    vec_dir_fy[jxi * nsph * _num_theta * _num_phi + jth * nsph * _num_phi + jph * nsph + i],
		    vec_dir_fz[jxi * nsph * _num_theta * _num_phi + jth * nsph * _num_phi + jph * nsph + i]
		  );
		}
		fprintf(p_outfile, "  MULS\n");
		for (int j = 0; j < 4; j++) {
		  int muls_index = 16 * cindex + 4 * j;
		  fprintf(
		    p_outfile,  "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
		    vec_dir_muls[muls_index],
		    vec_dir_muls[muls_index + 1],
		    vec_dir_muls[muls_index + 2],
		    vec_dir_muls[muls_index + 3]
		  );
		} // j muls loop
		fprintf(p_outfile, "  MULSLR\n");
		for (int j = 0; j < 4; j++) {
		  int muls_index = 16 * cindex + 4 * j;
		  fprintf(
		    p_outfile,  "        %15.7lE%15.7lE%15.7lE%15.7lE\n",
		    vec_dir_mulslr[muls_index],
		    vec_dir_mulslr[muls_index + 1],
		    vec_dir_mulslr[muls_index + 2],
		    vec_dir_mulslr[muls_index + 3]
		  );
		} // j mulslr loop
	      } // i sphere loop
	      done_dirs++;
	    } // jphs loop
	  } // jths loop
	} // jph loop
      } // jth loop
    } // jxi wavelength loop
    fclose(p_outfile);
  } else {
    result = -1;
  }
  return result;
}

#ifdef MPI_VERSION
int SphereOutputInfo::mpireceive(const mixMPI *mpidata, int pid) {
  int result = 0;
  int flag;
  int chk_nsph, chk_inpol, chk_isam, chk_num_theta, chk_num_thetas;
  int chk_num_phi, chk_num_phis, chk_ndirs, chk_idfc, chk_configs;
  double chk_exri;
  MPI_Recv(&flag, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Proceed with the rest _only if_ flag==0, else nothing is to be received
  if (flag == 0) {
    MPI_Recv(&chk_nsph, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_inpol, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_isam, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_theta, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_thetas, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_phi, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_num_phis, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_ndirs, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_exri, 1, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_idfc, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&chk_configs, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    result += (chk_nsph == nsph) ? 0 : 1;
    result += (chk_inpol == inpol) ? 0 : 1;
    result += (chk_isam == isam) ? 0 : 1;
    result += (chk_num_theta == _num_theta) ? 0 : 1;
    result += (chk_num_thetas == _num_thetas) ? 0 : 1;
    result += (chk_num_phi == _num_phi) ? 0 : 1;
    result += (chk_num_phis == _num_phis) ? 0 : 1;
    result += (chk_ndirs == ndirs) ? 0 : 1;
    result += (chk_exri == exri) ? 0 : 1;
    result += (chk_idfc == idfc) ? 0 : 1;
    result += (chk_configs == configurations) ? 0 : 1;
    if (result == 0) {
      int xi1, offset, chunk_size;
      MPI_Send(&result, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD);
      MPI_Recv(&xi1, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Receive vectors of single values per scale
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = xi1 - _first_xi;
      MPI_Recv(vec_jxi + offset, chunk_size, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_ier + offset, chunk_size, MPI_SHORT, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_vk + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_xi + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (nsph != 1) {
	MPI_Recv(vec_fsat + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
	MPI_Recv(vec_qschut + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(vec_pschut + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(vec_s0magt + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      // Receive vectors whose sizes depend on configurations and wavelengths
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * configurations;
      MPI_Recv(vec_sphere_sizes + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_sphere_ref_indices + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_abs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_albeds + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_scsrt + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_absrt + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_exsrt + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_fsas + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_qschu + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_pschu + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_s0mag + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_cosav + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_raprs + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqek1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqek2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqsk1 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_tqsk2 + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors whose sizes depend on NSPH, directions and wavelengths
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * nsph * ndirs;
      MPI_Recv(vec_dir_sas11 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas21 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas12 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_sas22 + offset, chunk_size, MPI_C_DOUBLE_COMPLEX, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_muls + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_mulslr + 16 * offset, 16 * chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive vectors whose sizes depend on NSPH, incidence directions and wavelengths
      MPI_Recv(&chunk_size, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset = (xi1 - _first_xi) * nsph * _num_theta * _num_phi;
      MPI_Recv(vec_dir_fx + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fy + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vec_dir_fz + offset, chunk_size, MPI_DOUBLE, pid, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }  else {
      MPI_Send(&result, 1, MPI_INT32_T, pid, 10, MPI_COMM_WORLD);
    }
  }
  return result;
}

int SphereOutputInfo::mpisend(const mixMPI *mpidata) {
  int result = 0;
  int chunk_size;
  if (_skip_flag == 1) {
    // tell the receiver we are not sending anything
    int flag = 1;
    MPI_Send(&flag, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
  }
  else {
    // tell the receiver we are sending actual stuff
    int flag = 0;
    MPI_Send(&flag, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    // Send output metadata for configuration cross-check
    MPI_Send(&nsph, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&inpol, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&isam, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_theta, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_thetas, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_phi, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&_num_phis, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&ndirs, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&exri, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&idfc, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    MPI_Send(&configurations, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
    // Wait for process 0 to cross-check the configuration
    MPI_Recv(&result, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (result == 0) {
      // Process 0 confirmed the consistency of configuration. Send the data.
      // Send vectors of single values per scale
      MPI_Send(&_first_xi, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(&xi_block_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_jxi, xi_block_size, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_ier, xi_block_size, MPI_SHORT, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_vk, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_xi, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      if (nsph != 1) {
	MPI_Send(vec_fsat, xi_block_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
	MPI_Send(vec_qschut, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
	MPI_Send(vec_pschut, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
	MPI_Send(vec_s0magt, xi_block_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      }

      // Send vectors whose sizes depend on configurations and scales
      chunk_size = xi_block_size * configurations;
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_sizes, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_sphere_ref_indices, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_abs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_albeds, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_scsrt, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_absrt, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_exsrt, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_fsas, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_qschu, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_pschu, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_s0mag, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_cosav, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_raprs, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqek1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqek2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqsk1, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_tqsk2, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);

      // Send vectors whose sizes depend on NSPH, directions and wavelengths
      chunk_size = xi_block_size * nsph * ndirs;
      MPI_Send(&chunk_size, 1, MPI_INT32_T, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas11, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas21, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas12, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_sas22, chunk_size, MPI_C_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_muls, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_mulslr, 16 * chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);

      // Send vectors whose sizes depend on NSPH, incidence directions and wavelengths
      chunk_size = xi_block_size * nsph * _num_theta * _num_phi;
      MPI_Send(vec_dir_fx, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fy, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
      MPI_Send(vec_dir_fz, chunk_size, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    }
  }
  return result;
}
#endif // MPI_VERSION
// >>> END OF SphereOutputInfo CLASS IMPLEMENTATION <<<
