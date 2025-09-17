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

/*! \file Commons.cpp
 *
 * \brief Implementation of the common data structures.
 */
#include <cstring>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H
#include "../include/Commons.h"
#endif

mixMPI::mixMPI() {
  mpirunning = 0;
  rank = 0;
  nprocs = 1;
}

#ifdef MPI_VERSION
mixMPI::mixMPI(MPI_Comm comm){
  mpirunning = 1;
  int ierr;
  // we should add some meaningful error checking and management here
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
}
#endif

mixMPI::mixMPI(const mixMPI& rhs) {
  mpirunning = rhs.mpirunning;
  rank = rhs.rank;
  nprocs = rhs.nprocs;
}

mixMPI::~mixMPI() {
}

// >>> ParticleDescriptor class implementation. <<< //
ParticleDescriptor::ParticleDescriptor(GeometryConfiguration *gconf, ScattererConfiguration *sconf) {
  _class_type = BASE_TYPE;
  _nsph = gconf->number_of_spheres;
  _li = (_nsph == 1) ? gconf->l_max : gconf->li;
  _le = 0;
  _lm = _li;
  _nlim = 0;
  _nlem = 0;
  _nlemt = 0;
  _ncou = 0;
  _litpo = 0;
  _litpos = 0;
  _lmpo = 0;
  _lmtpo = 0;
  _lmtpos = 0;
  _nv3j = 0;
  _ndi = 0;
  _ndit = 0;
  _ndm = 0;
  gcs = 0.0;
  _num_configurations = sconf->configurations;
  // _num_layers = (sconf->use_external_sphere) ? 1 : 0;
  _num_layers = 0;
  _max_layers = 1;
  for (int nli = 0; nli < num_configurations; nli++) {
    int nl = sconf->get_nshl(nli);
    if (nli == 0 && sconf->use_external_sphere) nl++;
    _num_layers += nl;
    if (nl >  _max_layers) _max_layers = nl;
  }

  vec_rmi = new dcomplex[_li * _nsph]();
  rmi = new dcomplex*[_li];
  vec_rei = new dcomplex[_li * _nsph]();
  rei = new dcomplex*[_li];
  for (int ri = 0; ri < _li; ri++) {
    rmi[ri] = vec_rmi + (_nsph * ri);
    rei[ri] = vec_rei + (_nsph * ri);
  }
  // int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  // vec_w = new dcomplex[nllt * 4]();
  // w = new dcomplex*[nllt];
  // for (int wi = 0; wi < nllt; wi++) w[wi] = vec_w + (4 * wi);
  vec_w = NULL;
  w = NULL;
  vec_rc = new double[num_layers]();
  rc = new double*[num_configurations];
  int last_layer_index = 0, cur_layer_index = 0;
  for (int rci = 0; rci < num_configurations; rci++) {
    rc[rci] = vec_rc + last_layer_index;
    last_layer_index += sconf->get_nshl(rci);
    if (rci == 0 && sconf->use_external_sphere) last_layer_index++;
    for (int rcj = 0; rcj < last_layer_index - cur_layer_index; rcj++) {
      rc[rci][rcj] = sconf->get_rcf(rci, rcj);
    }
    cur_layer_index = last_layer_index;
  }
  vint = new dcomplex[16]();
  rxx = new double[nsph]();
  ryy = new double[nsph]();
  rzz = new double[nsph]();
  iog = new int[nsph]();
  for (int si = 0; si < nsph; si++) {
    rxx[si] = gconf->get_sph_x(si);
    ryy[si] = gconf->get_sph_y(si);
    rzz[si] = gconf->get_sph_z(si);
    iog[si] = sconf->get_iog(si);
  }
  ros = new double[num_configurations]();
  nshl = new int[num_configurations]();
  for (int ci = 0; ci < num_configurations; ci++) {
    ros[ci] = sconf->get_radius(ci);
    nshl[ci] = sconf->get_nshl(ci);
  }
  _npnt = gconf->npnt;
  _npntts = gconf->npntts;
  int max_n = (npnt > npntts) ? npnt : npntts;
  _nhspo = 2 * max_n - 1;
  ris = new dcomplex[_nhspo]();
  dlri = new dcomplex[_nhspo]();
  vkt = new dcomplex[_nsph]();
  dc0 = new dcomplex[_max_layers + 1]();
  vsz = new double[_nsph]();
  
  // >>> NEEDED BY SPHERE AND CLUSTER <<<
  sas = NULL;
  vints = NULL;
  fsas = NULL;
  sscs = NULL;
  sexs = NULL;
  sabs = NULL;
  sqscs = NULL;
  sqexs = NULL;
  sqabs = NULL;
  gcsv = NULL;
  // >>> NEEDED BY CLUSTER <<<
  vec_tsas = NULL;
  vintt = NULL;
  tfsas = 0.0 + I * 0.0;
  tsas = NULL;
  gcs = 0.0;
  scs = 0.0;
  ecs = 0.0;
  acs = 0.0;
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  vec_am0m = NULL;
  vec_fsac = NULL;
  vec_sac = NULL;
  vec_fsacm = NULL;
  vec_ind3j = NULL;
  vh = NULL;
  vj0 = NULL;
  vyhj = NULL;
  vyj0 = NULL;
  vj = 0.0 + I * 0.0;
  am0m = NULL;
  fsac = NULL;
  sac = NULL;
  fsacm = NULL;
  vintm = NULL;
  scscp = NULL;
  ecscp = NULL;
  scscpm = NULL;
  ecscpm = NULL;
  v3j0 = NULL;
  scsc = NULL;
  ecsc = NULL;
  scscm = NULL;
  ecscm = NULL;
  ind3j = NULL;
  rac3j = NULL;
  // >>> NEEDED BY INCLU <<<
  rm0 = NULL;
  re0 = NULL;
  rmw = NULL;
  rew = NULL;
  tm = NULL;
  te = NULL;
  tm0 = NULL;
  te0 = NULL;
  vec_at = NULL;
  at = NULL;
}

ParticleDescriptor::ParticleDescriptor(const ParticleDescriptor &rhs) {
  _class_type = rhs._class_type;
  _nsph = rhs._nsph;
  _li = rhs._li;
  _le = rhs._le;
  _lm = rhs._lm;
  _nlim = rhs._nlim;
  _nlem = rhs._nlem;
  _nlemt = rhs._nlemt;
  _ncou = rhs._ncou;
  _litpo = rhs._litpo;
  _litpos = rhs._litpos;
  _lmpo = rhs._lmpo;
  _lmtpo = rhs._lmtpo;
  _lmtpos = rhs._lmtpos;
  _nv3j = rhs._nv3j;
  _ndi = rhs._ndi;
  _ndit = rhs._ndit;
  _ndm = rhs._ndm;
  gcs = rhs.gcs;
  _max_layers = rhs._max_layers;
  _num_configurations = rhs._num_configurations;
  _num_layers = rhs._num_layers;

  vec_rmi = new dcomplex[_li * _nsph];
  rmi = new dcomplex*[_li];
  vec_rei = new dcomplex[_li * _nsph];
  rei = new dcomplex*[_li];
  for (int ri = 0; ri < _li; ri++) {
    rmi[ri] = vec_rmi + (_nsph * ri);
    rei[ri] = vec_rei + (_nsph * ri);
  }
  for (int rj = 0; rj < _li * _nsph; rj++) {
    vec_rmi[rj] = rhs.vec_rmi[rj];
    vec_rei[rj] = rhs.vec_rei[rj];
  }
  int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  vec_w = new dcomplex[nllt * 4];
  w = new dcomplex*[nllt];
  for (int wi = 0; wi < nllt; wi++) w[wi] = vec_w + (4 * wi);
  for (int wj = 0; wj < 4 * nllt; wj++) vec_w[wj] = rhs.vec_w[wj];
  vec_rc = new double[_num_layers];
  for (int rcj = 0; rcj < _num_layers; rcj++) vec_rc[rcj] = rhs.vec_rc[rcj];
  int tot_nshl = 0;
  for (int nsi = 0; nsi < _num_configurations; nsi++) tot_nshl += rhs.nshl[nsi];
  rc = new double*[_num_configurations];
  int last_layer_index = 0;
  for (int rci = 0; rci < _num_configurations; rci++) {
    rc[rci] = vec_rc + last_layer_index;
    last_layer_index += rhs.nshl[rci];
    if (rci == 0 && tot_nshl < _num_layers) last_layer_index++;
  }
  vint = new dcomplex[16];
  for (int vi = 0; vi < 16; vi++) vint[vi] = rhs.vint[vi];
  rxx = new double[_nsph];
  ryy = new double[_nsph];
  rzz = new double[_nsph];
  iog = new int[_nsph];
  for (int si = 0; si < _nsph; si++) {
    rxx[si] = rhs.rxx[si];
    ryy[si] = rhs.ryy[si];
    rzz[si] = rhs.rzz[si];
    iog[si] = rhs.iog[si];
  }
  ros = new double[_num_configurations];
  nshl = new int[_num_configurations];
  for (int ci = 0; ci < _num_configurations; ci++) {
    ros[ci] = rhs.ros[ci];
    nshl[ci] = rhs.nshl[ci];
  }
  _npnt = rhs._npnt;
  _npntts = rhs._npntts;
  _nhspo = rhs._nhspo;
  _max_layers = rhs._max_layers;
  ris = new dcomplex[_nhspo]();
  dlri = new dcomplex[_nhspo]();
  for (int ri = 0; ri < _nhspo; ri++) {
    ris[ri] = rhs.ris[ri];
    dlri[ri] = rhs.dlri[ri];
  }
  vkt = new dcomplex[_nsph]();
  vsz = new double[_nsph]();
  for (int vi = 0; vi < _nsph; vi++) {
    vkt[vi] = rhs.vkt[vi];
    vsz[vi] = rhs.vsz[vi];
  }
  dc0 = new dcomplex[_max_layers]();
  for (int di = 0; di < _max_layers; di++) {
    dc0[di] = rhs.dc0[di];
  }
  // >>> NEEDED BY SPHERE AND CLUSTER <<<
  sas = NULL;
  vints = NULL;
  fsas = NULL;
  sscs = NULL;
  sexs = NULL;
  sabs = NULL;
  sqscs = NULL;
  sqexs = NULL;
  sqabs = NULL;
  gcsv = NULL;
  // >>> NEEDED BY CLUSTER <<<
  vintt = NULL;
  tfsas = 0.0 + I * 0.0;
  vec_tsas = NULL;
  tsas = NULL;
  scs = 0.0;
  ecs = 0.0;
  acs = 0.0;
  vec_gis = NULL;
  vec_gls = NULL;
  vec_sam = NULL;
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  vec_am0m = NULL;
  vec_fsac = NULL;
  vec_sac = NULL;
  vec_fsacm = NULL;
  vec_ind3j = NULL;
  vh = NULL;
  vj0 = NULL;
  vyhj = NULL;
  vyj0 = NULL;
  vj = 0.0 + I * 0.0;
  am0m = NULL;
  fsac = NULL;
  sac = NULL;
  fsacm = NULL;
  vintm = NULL;
  scscp = NULL;
  ecscp = NULL;
  scscpm = NULL;
  ecscpm = NULL;
  v3j0 = NULL;
  scsc = NULL;
  ecsc = NULL;
  scscm = NULL;
  ecscm = NULL;
  ind3j = NULL;
  rac3j = NULL;
  // >>> NEEDED BY INCLU <<<
  rm0 = NULL;
  re0 = NULL;
  rmw = NULL;
  rew = NULL;
  tm = NULL;
  te = NULL;
  tm0 = NULL;
  te0 = NULL;
}

#ifdef MPI_VERSION
ParticleDescriptor::ParticleDescriptor(const mixMPI *mpidata) {
  MPI_Bcast(&_class_type, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nsph, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_li, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_le, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nlim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nlem, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nlemt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ncou, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_litpo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_litpos, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lmpo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lmtpo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lmtpos, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nv3j, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ndi, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ndit, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ndm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gcs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_num_configurations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_num_layers, 1, MPI_INT, 0, MPI_COMM_WORLD);

  vec_rmi = new dcomplex[_li * _nsph];
  vec_rei = new dcomplex[_li * _nsph];
  MPI_Bcast(vec_rmi, _nsph * _li, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_rei, _nsph * _li, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  rmi = new dcomplex*[_li];
  rei = new dcomplex*[_li];
  for (int ri = 0; ri < _li; ri++) {
    rmi[ri] = vec_rmi + (_nsph * ri);
    rei[ri] = vec_rei + (_nsph * ri);
  }
  int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  vec_w = new dcomplex[nllt * 4];
  MPI_Bcast(vec_w, 4 * nllt, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  w = new dcomplex*[nllt];
  for (int wi = 0; wi < nllt; wi++) w[wi] = vec_w + (4 * wi);
  vec_rc = new double[_num_layers];
  MPI_Bcast(vec_rc, _num_layers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  nshl = new int[num_configurations];
  MPI_Bcast(nshl, _num_configurations, MPI_INT, 0, MPI_COMM_WORLD);
  int tot_nshl = 0;
  for (int nsi = 0; nsi < _num_configurations; nsi++) tot_nshl += nshl[nsi];
  rc = new double*[_num_configurations];
  int last_layer_index = 0;
  for (int rci = 0; rci < _num_configurations; rci++) {
    rc[rci] = vec_rc + last_layer_index;
    last_layer_index += nshl[rci];
    if (rci == 0 && tot_nshl < _num_layers) last_layer_index++;
  }
  vint = new dcomplex[16];
  MPI_Bcast(vint, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  rxx = new double[_nsph];
  ryy = new double[_nsph];
  rzz = new double[_nsph];
  iog = new int[_nsph];
  MPI_Bcast(rxx, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ryy, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(rzz, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(iog, _nsph, MPI_INT, 0, MPI_COMM_WORLD);
  ros = new double[_num_configurations];
  MPI_Bcast(ros, _num_configurations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npntts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nhspo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_max_layers, 1, MPI_INT, 0, MPI_COMM_WORLD);
  ris = new dcomplex[_nhspo];
  MPI_Bcast(ris, _nhspo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  dlri = new dcomplex[_nhspo];
  MPI_Bcast(dlri, _nhspo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vkt = new dcomplex[_nsph];
  MPI_Bcast(vkt, _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vsz = new double[_nsph];
  MPI_Bcast(vsz, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  dc0 = new dcomplex[_max_layers];
  MPI_Bcast(dc0, _max_layers, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}

void ParticleDescriptor::mpibcast(const mixMPI *mpidata) {
  MPI_Bcast(&_class_type, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nsph, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_li, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_le, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nlim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nlem, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nlemt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ncou, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_litpo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_litpos, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lmpo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lmtpo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_lmtpos, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nv3j, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ndi, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ndit, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ndm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gcs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_num_configurations, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_num_layers, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_rmi, _nsph * _li, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_rei, _nsph * _li, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  MPI_Bcast(vec_w, 4 * nllt, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vec_rc, _num_layers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(nshl, _num_configurations, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vint, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(rxx, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ryy, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(rzz, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(iog, _nsph, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ros, _num_configurations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npnt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_npntts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nhspo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_max_layers, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ris, _nhspo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(dlri, _nhspo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vkt, _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(vsz, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(dc0, _max_layers, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  // >>> NEEDED BY SPHERE AND CLUSTER <<< //
  if (_class_type == SPHERE_TYPE || _class_type == CLUSTER_TYPE) {
    MPI_Bcast(vec_sas, 4 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_vints, 16 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(fsas, _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(sscs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sexs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sabs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sqscs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sqexs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sqabs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(gcsv, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  // >>> NEEDED BY CLUSTER <<< //
  if (_class_type == CLUSTER_TYPE) {
    MPI_Bcast(vintt, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_tsas, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tfsas, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&scs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ecs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&acs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_gis, _ndi * _nlem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_gls, _ndi * _nlem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_sam, _ndit * _nlemt, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  }
  // >>> NEEDED BY CLUSTER AND INCLU <<< //
  if (_class_type == CLUSTER_TYPE || _class_type == INCLUSION_TYPE) {
    int nlemts = _nlemt * _nlemt;
    MPI_Bcast(vec_am0m, nlemts, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_fsac, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_sac, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_fsacm, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    int vec_ind3j_size = (_lm + 1) * _lm;
    MPI_Bcast(vec_ind3j, vec_ind3j_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(vh, _ncou * _litpo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vj0, _nsph * _lmtpo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vyhj, _ncou * _litpos, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vyj0, _nsph * _lmtpos, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vintm, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(scscp, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(ecscp, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(scscpm, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(ecscpm, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(scsc, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ecsc, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(scscm, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ecscm, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v3j0, _nv3j, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rac3j, _lmtpo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  // >>> NEEDED BY INCLU <<< //
  if (_class_type == INCLUSION_TYPE) {
    MPI_Bcast(rm0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(re0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(rmw, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(rew, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tm, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(te, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(tm0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(te0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(vec_at, _nlemt * _ndm, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  }
}
#endif // MPI_VERSION

ParticleDescriptor::~ParticleDescriptor() {
  // Base class members, always destroyed.
  delete[] nshl;
  delete[] ros;
  delete[] iog;
  delete[] rzz;
  delete[] ryy;
  delete[] rxx;
  delete[] vint;
  delete[] rc;
  delete[] vec_rc;
  delete[] w;
  delete[] vec_w;
  delete[] rei;
  delete[] vec_rei;
  delete[] rmi;
  delete[] vec_rmi;
  delete[] ris;
  delete[] dlri;
  delete[] vkt;
  delete[] dc0;
  delete[] vsz;
  // Inclusion class members, destroyed only if sub-class is INCLUSION
  if (_class_type == INCLUSION_TYPE) {
    delete[] rm0;
    delete[] re0;
    delete[] rmw;
    delete[] rew;
    delete[] tm;
    delete[] te;
    delete[] tm0;
    delete[] te0;
    delete[] at;
    delete[] vec_at;
  }
  // Inclusion/cluster class members, destroyed if sub-class is INCLUSION or CLUSTER
  if (_class_type == INCLUSION_TYPE || _class_type == CLUSTER_TYPE) {
    delete[] vec_am0m;
    delete[] vec_fsac;
    delete[] vec_sac;
    delete[] vec_fsacm;
    delete[] vec_ind3j;
  
    delete[] vh;
    delete[] vj0;
    delete[] vyhj;
    delete[] vyj0;
    delete[] am0m;
    delete[] sac;
    delete[] fsac;
    delete[] fsacm;
    delete[] vintm;
    delete[] scscp;
    delete[] ecscp;
    delete[] scscpm;
    delete[] ecscpm;
    delete[] v3j0;
    delete[] scsc;
    delete[] ecsc;
    delete[] scscm;
    delete[] ecscm;
    delete[] ind3j;
    delete[] rac3j;
  }
  // Cluster/sphere class members, destroyed if sub-class is CLUSTER or SPHERE
  if (_class_type == CLUSTER_TYPE || _class_type == SPHERE_TYPE) {
    for (int vi = 0; vi < _nsph; vi++) {
      delete[] sas[vi];
    }
    delete[] vints;
    delete[] sas;
    delete[] gcsv;
    delete[] sqabs;
    delete[] sqexs;
    delete[] sqscs;
    delete[] sabs;
    delete[] sexs;
    delete[] sscs;
    delete[] fsas;
    delete[] vec_vints;
    delete[] vec_sas;
  }
  // Cluster class members, destroyed only if sub-class is CLUSTER
  if (_class_type == CLUSTER_TYPE) {
    delete[] vintt;
    delete[] vec_tsas;
    delete[] tsas;
    delete[] vec_gis;
    delete[] vec_gls;
    delete[] vec_sam;
    delete[] gis;
    delete[] gls;
    delete[] sam;
  }
}

// >>> End of ParticleDescriptor class implementation. <<< //

// >>> ParticleDescriptorCluster class implementation. <<< //
ParticleDescriptorCluster::ParticleDescriptorCluster(GeometryConfiguration *gconf, ScattererConfiguration *sconf) : ParticleDescriptor(gconf, sconf) {
  _class_type = CLUSTER_TYPE;
  // Needed by SPHERE and CLUSTER
  vec_sas = new dcomplex[4 * _nsph]();
  vec_vints = new dcomplex[16 * _nsph]();

  fsas = new dcomplex[_nsph];
  sscs = new double[_nsph]();
  sexs = new double[_nsph]();
  sabs = new double[_nsph]();
  sqscs = new double[_nsph]();
  sqexs = new double[_nsph]();
  sqabs = new double[_nsph]();
  gcsv = new double[_nsph]();
  sas = new dcomplex**[_nsph];
  vints = new dcomplex*[_nsph];
  for (int vi = 0; vi < _nsph; vi++) {
    vints[vi] = vec_vints + (16 * vi);
    sas[vi] = new dcomplex*[2];
    sas[vi][0] = vec_sas + (4 * vi);
    sas[vi][1] = vec_sas + (4 * vi) + 2;
  }

  // Needed by CLUSTER
  vintt = new dcomplex[16]();
  vec_tsas = new dcomplex[4]();
  tsas = new dcomplex*[2];
  tsas[0] = vec_tsas;
  tsas[1] = vec_tsas + 2;

  // Needed by CLUSTER and INCLU
  _le = gconf->le;
  _lm = (_li > _le) ? _li : _le;
  _nlim = _li * (_li + 2);
  _nlem = _le * (_le + 2);
  _nlemt = 2 * _nlem;
  _ncou = _nsph * _nsph - 1;
  _litpo = _li + _li + 1;
  _litpos = _litpo * _litpo;
  _lmpo = _lm + 1;
  _lmtpo = _li + _le + 1;
  _lmtpos = _lmtpo * _lmtpo;
  _nv3j = (_lm * (_lm + 1) * (2 * _lm + 7)) / 6;
  _ndi = _nsph * _nlim;
  _ndit = 2 * _nsph * _nlim;
  int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  vec_w = new dcomplex[nllt * 4]();
  w = new dcomplex*[nllt];
  for (int wi = 0; wi < nllt; wi++) w[wi] = vec_w + (4 * wi);
  vec_am0m = new dcomplex[_nlemt * _nlemt]();
  vec_fsac = new dcomplex[4]();
  vec_sac = new dcomplex[4]();
  vec_fsacm = new dcomplex[4]();
  vec_ind3j = new int[(_lm + 1) * _lm]();
  
  vh = new dcomplex[_ncou * _litpo]();
  vj0 = new dcomplex[_nsph * _lmtpo]();
  vyhj = new dcomplex[_ncou * _litpos]();
  vyj0 = new dcomplex[_nsph * _lmtpos]();
  am0m = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
  sac = new dcomplex*[2];
  fsac = new dcomplex*[2];
  fsacm = new dcomplex*[2];
  for (int fi = 0; fi < 2; fi++) {
    sac[fi] = vec_sac + (fi * 2);
    fsac[fi] = vec_fsac + (fi * 2);
    fsacm[fi] = vec_fsacm + (fi * 2);
  }
  vintm = new dcomplex[16]();
  scscp = new dcomplex[2]();
  ecscp = new dcomplex[2]();
  scscpm = new dcomplex[2]();
  ecscpm = new dcomplex[2]();
  v3j0 = new double[_nv3j]();
  scsc = new double[2]();
  ecsc = new double[2]();
  scscm = new double[2]();
  ecscm = new double[2]();
  ind3j = new int*[_lm + 1];
  for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
  rac3j = new double[_lmtpo]();

  // Needed by CLUSTER
  vec_gis = new dcomplex[_ndi * _nlem]();
  gis = new dcomplex*[_ndi];
  vec_gls = new dcomplex[_ndi * _nlem]();
  gls = new dcomplex*[_ndi];
  for (int gi = 0; gi < _ndi; gi++) {
    gis[gi] = vec_gis + (gi * _nlem);
    gls[gi] = vec_gls + (gi * _nlem);
  }
  vec_sam = new dcomplex[_ndit * _nlemt]();
  sam = new dcomplex*[_ndit];
  for (int si = 0; si < _ndit; si++) sam[si] = vec_sam + (si * _nlemt);
}

ParticleDescriptorCluster::ParticleDescriptorCluster(const ParticleDescriptorCluster &rhs) : ParticleDescriptor(rhs) {
  // >>> NEEDED BY SPHERE AND CLUSTER <<<
  vec_sas = new dcomplex[4 * _nsph];
  for (int vsi = 0; vsi < 4 * _nsph; vsi++) vec_sas[vsi] = rhs.vec_sas[vsi];
  vec_vints = new dcomplex[16 * _nsph];
  for (int vvi = 0; vvi < 16 * _nsph; vvi++) vec_vints[vvi] = rhs.vec_vints[vvi];
  
  fsas = new dcomplex[_nsph];
  sscs = new double[_nsph];
  sexs = new double[_nsph];
  sabs = new double[_nsph];
  sqscs = new double[_nsph];
  sqexs = new double[_nsph];
  sqabs = new double[_nsph];
  gcsv = new double[_nsph];
  for (int gi = 0; gi < _nsph; gi++) {
    fsas[gi] = rhs.fsas[gi];
    sscs[gi] = rhs.sscs[gi];
    sexs[gi] = rhs.sexs[gi];
    sabs[gi] = rhs.sabs[gi];
    sqscs[gi] = rhs.sqscs[gi];
    sqexs[gi] = rhs.sqexs[gi];
    sqabs[gi] = rhs.sqabs[gi];
    gcsv[gi] = rhs.gcsv[gi];
  }
  sas = new dcomplex**[_nsph];
  vints = new dcomplex*[_nsph];
  for (int vi = 0; vi < nsph; vi++) {
    vints[vi] = vec_vints + (16 * vi);
    sas[vi] = new dcomplex*[2];
    sas[vi][0] = vec_sas + (4 * vi);
    sas[vi][1] = vec_sas + (4 * vi) + 2;
  }
  // >>> NEEDED BY CLUSTER <<<
  vintt = new dcomplex[16];
  for (int ti = 0; ti < 16; ti++) vintt[ti] = rhs.vintt[ti];
  vec_tsas = new dcomplex[4];
  for (int si = 0; si < 4; si++) vec_tsas[si] = rhs.vec_tsas[si];
  tsas = new dcomplex*[2];
  tsas[0] = vec_tsas;
  tsas[1] = vec_tsas + 2;
  tfsas = rhs.tfsas;
  gcs = rhs.gcs;
  scs = rhs.scs;
  ecs = rhs.ecs;
  acs = rhs.acs;
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  vec_am0m = new dcomplex[_nlemt * _nlemt];
  np_int nlemts = _nlemt * _nlemt;
  for (np_int ai = 0; ai < nlemts; ai++) vec_am0m[ai] = rhs.vec_am0m[ai];
  vec_fsac = new dcomplex[4];
  vec_sac = new dcomplex[4];
  vec_fsacm = new dcomplex[4];
  for (int vj = 0; vj < 4; vj++) {
    vec_fsac[vj] = rhs.vec_fsac[vj];
    vec_sac[vj] = rhs.vec_sac[vj];
    vec_fsacm[vj] = rhs.vec_fsacm[vj];
  }
  vec_ind3j = new int[(_lm + 1) * _lm];
  np_int vec_ind3j_size = (_lm + 1) * _lm;
  for (np_int ii = 0; ii < vec_ind3j_size; ii++) vec_ind3j[ii] = rhs.vec_ind3j[ii];
  
  vh = new dcomplex[_ncou * _litpo];
  for (int hi = 0; hi < _ncou * _litpo; hi++) vh[hi] = rhs.vh[hi];
  vj0 = new dcomplex[_nsph * _lmtpo];
  for (int ji = 0; ji < _nsph * _lmtpo; ji++) vj0[ji] = rhs.vj0[ji];
  vyhj = new dcomplex[_ncou * _litpos];
  for (int hj = 0; hj < _ncou * _litpos; hj++) vyhj[hj] = rhs.vyhj[hj];
  vyj0 = new dcomplex[_nsph * _lmtpos];
  for (int jj = 0; jj < _nsph * _lmtpos; jj++) vyj0[jj] = rhs.vyj0[jj];
  am0m = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
  sac = new dcomplex*[2];
  fsac = new dcomplex*[2];
  fsacm = new dcomplex*[2];
  for (int fi = 0; fi < 2; fi++) {
    sac[fi] = vec_sac + (fi * 2);
    fsac[fi] = vec_fsac + (fi * 2);
    fsacm[fi] = vec_fsacm + (fi * 2);
  }
  vintm = new dcomplex[16];
  for (int mi = 0; mi < 16; mi++) vintm[mi] = rhs.vintm[mi];
  scscp = new dcomplex[2];
  ecscp = new dcomplex[2];
  scscpm = new dcomplex[2];
  ecscpm = new dcomplex[2];
  scsc = new double[2];
  ecsc = new double[2];
  scscm = new double[2];
  ecscm = new double[2];
  for (int ci = 0; ci < 2; ci++) {
    scscp[ci] = rhs.scscp[ci];
    ecscp[ci] = rhs.ecscp[ci];
    scscpm[ci] = rhs.scscpm[ci];
    ecscpm[ci] = rhs.ecscpm[ci];
    scsc[ci] = rhs.scsc[ci];
    ecsc[ci] = rhs.ecsc[ci];
    scscm[ci] = rhs.scscm[ci];
    ecscm[ci] = rhs.ecscm[ci];
  }
  v3j0 = new double[_nv3j];
  // for (int vj = 0; vj < _nv3j; vj++) v3j0[vj] = rhs.v3j0[_nv3j]; REPORT: AAAAH! ORRORE E DISGUSTO!
  for (int vj = 0; vj < _nv3j; vj++) v3j0[vj] = rhs.v3j0[vj];
  ind3j = new int*[_lm + 1];
  for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
  rac3j = new double[_lmtpo];
  for (int ri = 0; ri < _lmtpo; ri++) rac3j[ri] = rhs.rac3j[ri];

  // Needed by CLUSTER
  vec_gis = new dcomplex[_ndi * _nlem];
  vec_gls = new dcomplex[_ndi * _nlem];
  for (int vi = 0; vi < _ndi * _nlem; vi++) {
    vec_gis[vi] = rhs.vec_gis[vi];
    vec_gls[vi] = rhs.vec_gls[vi];
  }
  gis = new dcomplex*[_ndi];
  gls = new dcomplex*[_ndi];
  for (int gi = 0; gi < _ndi; gi++) {
    gis[gi] = vec_gis + (gi * _nlem);
    gls[gi] = vec_gls + (gi * _nlem);
  }
  vec_sam = new dcomplex[_ndit * _nlemt];
  for (int vi = 0; vi < _ndit * _nlemt; vi++) vec_sam[vi] = rhs.vec_sam[vi];
  sam = new dcomplex*[_ndit];
  for (int si = 0; si < _ndit; si++) sam[si] = vec_sam + (si * _nlemt);
}

#ifdef MPI_VERSION
ParticleDescriptorCluster::ParticleDescriptorCluster(const mixMPI *mpidata) : ParticleDescriptor(mpidata) {
  // >>> NEEDED BY SPHERE AND CLUSTER <<< //
  vec_sas = new dcomplex[4 * _nsph];
  MPI_Bcast(vec_sas, 4 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_vints = new dcomplex[16 * _nsph];
  MPI_Bcast(vec_vints, 16 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  fsas = new dcomplex[_nsph];
  MPI_Bcast(fsas, _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  sscs = new double[_nsph];
  MPI_Bcast(sscs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sexs = new double[_nsph];
  MPI_Bcast(sexs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sabs = new double[_nsph];
  MPI_Bcast(sabs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sqscs = new double[_nsph];
  MPI_Bcast(sqscs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sqexs = new double[_nsph];
  MPI_Bcast(sqexs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sqabs = new double[_nsph];
  MPI_Bcast(sqabs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  gcsv = new double[_nsph];
  MPI_Bcast(gcsv, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sas = new dcomplex**[_nsph];
  vints = new dcomplex*[_nsph];
  for (int vi = 0; vi < nsph; vi++) {
    vints[vi] = vec_vints + (16 * vi);
    sas[vi] = new dcomplex*[2];
    sas[vi][0] = vec_sas + (4 * vi);
    sas[vi][1] = vec_sas + (4 * vi) + 2;
  }
  // >>> NEEDED BY CLUSTER <<< //
  vintt = new dcomplex[16];
  MPI_Bcast(vintt, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_tsas = new dcomplex[4];
  MPI_Bcast(vec_tsas, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  tsas = new dcomplex*[2];
  tsas[0] = vec_tsas;
  tsas[1] = vec_tsas + 2;
  MPI_Bcast(&tfsas, 1, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  MPI_Bcast(&scs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ecs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&acs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  vec_gis = new dcomplex[_ndi * _nlem];
  MPI_Bcast(vec_gis, _ndi * _nlem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_gls = new dcomplex[_ndi * _nlem];
  MPI_Bcast(vec_gls, _ndi * _nlem, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_sam = new dcomplex[_ndit * _nlemt];
  MPI_Bcast(vec_sam, _ndit * _nlemt, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  gis = new dcomplex*[_ndi];
  gls = new dcomplex*[_ndi];
  for (int gi = 0; gi < _ndi; gi++) {
    gis[gi] = vec_gis + (gi * _nlem);
    gls[gi] = vec_gls + (gi * _nlem);
  }
  sam = new dcomplex*[_ndit];
  for (int si = 0; si < _ndit; si++) {
    sam[si] = vec_sam + (si * _nlemt);
  }
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  vec_am0m = new dcomplex[_nlemt * _nlemt];
  int nlemts = _nlemt * _nlemt;
  MPI_Bcast(vec_am0m, nlemts, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_fsac = new dcomplex[4];
  MPI_Bcast(vec_fsac, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_sac = new dcomplex[4];
  MPI_Bcast(vec_sac, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_fsacm = new dcomplex[4];
  MPI_Bcast(vec_fsacm, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  int vec_ind3j_size = (_lm + 1) * _lm;
  vec_ind3j = new int[vec_ind3j_size];
  MPI_Bcast(vec_ind3j, vec_ind3j_size, MPI_INT, 0, MPI_COMM_WORLD);
  
  vh = new dcomplex[_ncou * _litpo];
  MPI_Bcast(vh, _ncou * _litpo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vj0 = new dcomplex[_nsph * _lmtpo];
  MPI_Bcast(vj0, _nsph * _lmtpo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vyhj = new dcomplex[_ncou * _litpos];
  MPI_Bcast(vyhj, _ncou * _litpos, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vyj0 = new dcomplex[_nsph * _lmtpos];
  MPI_Bcast(vyj0, _nsph * _lmtpos, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  am0m = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
  sac = new dcomplex*[2];
  fsac = new dcomplex*[2];
  fsacm = new dcomplex*[2];
  for (int fi = 0; fi < 2; fi++) {
    sac[fi] = vec_sac + (fi * 2);
    fsac[fi] = vec_fsac + (fi * 2);
    fsacm[fi] = vec_fsacm + (fi * 2);
  }
  vintm = new dcomplex[16];
  MPI_Bcast(vintm, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  scscp = new dcomplex[2];
  MPI_Bcast(scscp, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  ecscp = new dcomplex[2];
  MPI_Bcast(ecscp, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  scscpm = new dcomplex[2];
  MPI_Bcast(scscpm, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  ecscpm = new dcomplex[2];
  MPI_Bcast(ecscpm, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  scsc = new double[2];
  MPI_Bcast(scsc, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ecsc = new double[2];
  MPI_Bcast(ecsc, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  scscm = new double[2];
  MPI_Bcast(scscm, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ecscm = new double[2];
  MPI_Bcast(ecscm, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  v3j0 = new double[_nv3j];
  MPI_Bcast(v3j0, _nv3j, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ind3j = new int*[_lm + 1];
  for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
  rac3j = new double[_lmtpo];
  MPI_Bcast(rac3j, _lmtpo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
#endif // MPI_VERSION

int ParticleDescriptorCluster::update_orders(int inner_order, int outer_order) {
  int result = 0;
  bool changed_li = false;
  bool changed_le = false;
  if (inner_order != _li) {
    _li = inner_order;
    changed_li = true;
    delete[] vec_rmi;
    delete[] rmi;
    vec_rmi = new dcomplex[_li * _nsph]();
    rmi = new dcomplex*[_li];
    delete[] vec_rei;
    delete[] rei;
    vec_rei = new dcomplex[_li * _nsph]();
    rei = new dcomplex*[_li];
    for (int ri = 0; ri < _li; ri++) {
      rmi[ri] = vec_rmi + (_nsph * ri);
      rei[ri] = vec_rei + (_nsph * ri);
    }
    _litpo = _li + _li + 1;
    _litpos = _litpo * _litpo;
    delete[] vh;
    vh = new dcomplex[_ncou * _litpo]();
    delete[] vyhj;
    vyhj = new dcomplex[_ncou * _litpos]();
    _nlim = _li * (_li + 2);
    _ndi = _nsph * _nlim;
    _ndit = 2 * _nsph * _nlim;
  }
  if (outer_order != _le) {
    _le = outer_order;
    changed_le = true;
    _nlem = _le * (_le + 2);
    _nlemt = 2 * _nlem;
    delete[] vec_am0m;
    vec_am0m = new dcomplex[_nlemt * _nlemt]();
    delete[] am0m;
    am0m = new dcomplex*[_nlemt];
    for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
    delete[] vec_w;
    delete[] w;
    vec_w = new dcomplex[_nlemt * 4]();
    w = new dcomplex*[nlemt];
    for (int wi = 0; wi < nlemt; wi++) w[wi] = vec_w + (4 * wi);
  }
  if (changed_li || changed_le) {
    _lm = (_li > _le) ? _li : _le;
    _lmpo = _lm + 1;
    _lmtpo = _li + _le + 1;
    _lmtpos = _lmtpo * _lmtpo;
    _nv3j = (_lm * (_lm + 1) * (2 * _lm + 7)) / 6;
    delete[] vec_ind3j;
    vec_ind3j = new int[(_lm + 1) * _lm]();
    delete[] ind3j;
    ind3j = new int*[_lm + 1];
    for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
    delete[] vj0;
    vj0 = new dcomplex[_nsph * _lmtpo]();
    delete[] vyj0;
    vyj0 = new dcomplex[_nsph * _lmtpos]();
    delete[] v3j0;
    v3j0 = new double[_nv3j]();
    delete[] rac3j;
    rac3j = new double[_lmtpo]();
    delete[] vec_gis;
    vec_gis = new dcomplex[_ndi * _nlem]();
    delete[] gis;
    gis = new dcomplex*[_ndi];
    delete[] vec_gls;
    vec_gls = new dcomplex[_ndi * _nlem]();
    delete[] gls;
    gls = new dcomplex*[_ndi];
    for (int gi = 0; gi < _ndi; gi++) {
      gis[gi] = vec_gis + (gi * _nlem);
      gls[gi] = vec_gls + (gi * _nlem);
    }
    delete[] vec_sam;
    vec_sam = new dcomplex[_ndit * _nlemt]();
    delete[] sam;
    sam = new dcomplex*[_ndit];
    for (int si = 0; si < _ndit; si++) sam[si] = vec_sam + (si * _nlemt);
  }
  return result;
}
// >>> End of ParticleDescriptorCluster class implementation. <<< //

// >>> ParticleDescriptorInclusion class implementation. <<< //
ParticleDescriptorInclusion::ParticleDescriptorInclusion(GeometryConfiguration *gconf, ScattererConfiguration *sconf) : ParticleDescriptor(gconf, sconf) {
  _class_type = INCLUSION_TYPE;
  // Needed by CLUSTER and INCLU
  _le = gconf->le;
  _lm = (_li > _le) ? _li : _le;
  _nlim = _li * (_li + 2);
  _nlem = _le * (_le + 2);
  _nlemt = 2 * _nlem;
  _ncou = _nsph * _nsph - 1;
  _litpo = _li + _li + 1;
  _litpos = _litpo * _litpo;
  _lmpo = _lm + 1;
  _lmtpo = _li + _le + 1;
  _lmtpos = _lmtpo * _lmtpo;
  _nv3j = (_lm * (_lm + 1) * (2 * _lm + 7)) / 6;
  _ndi = _nsph * _nlim;
  _ndit = 2 * _nsph * _nlim;
  int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  vec_w = new dcomplex[nllt * 4]();
  w = new dcomplex*[nllt];
  for (int wi = 0; wi < nllt; wi++) w[wi] = vec_w + (4 * wi);
  vec_am0m = new dcomplex[_nlemt * _nlemt]();
  vec_fsac = new dcomplex[4]();
  vec_sac = new dcomplex[4]();
  vec_fsacm = new dcomplex[4]();
  vec_ind3j = new int[(_lm + 1) * _lm]();
  
  vh = new dcomplex[_ncou * _litpo]();
  vj0 = new dcomplex[_nsph * _lmtpo]();
  vyhj = new dcomplex[_ncou * _litpos]();
  vyj0 = new dcomplex[_nsph * _lmtpos]();
  am0m = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
  sac = new dcomplex*[2];
  fsac = new dcomplex*[2];
  fsacm = new dcomplex*[2];
  for (int fi = 0; fi < 2; fi++) {
    sac[fi] = vec_sac + (fi * 2);
    fsac[fi] = vec_fsac + (fi * 2);
    fsacm[fi] = vec_fsacm + (fi * 2);
  }
  vintm = new dcomplex[16]();
  scscp = new dcomplex[2]();
  ecscp = new dcomplex[2]();
  scscpm = new dcomplex[2]();
  ecscpm = new dcomplex[2]();
  v3j0 = new double[_nv3j]();
  scsc = new double[2]();
  ecsc = new double[2]();
  scscm = new double[2]();
  ecscm = new double[2]();
  ind3j = new int*[_lm + 1];
  for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
  rac3j = new double[_lmtpo]();
  // Needed by INCLU
  _ndm = 2 * (_nsph  * _nlim + _nlem);
  rm0 = new dcomplex[_le]();
  re0 = new dcomplex[_le]();
  rmw = new dcomplex[_le]();
  rew = new dcomplex[_le]();
  tm = new dcomplex[_le]();
  te = new dcomplex[_le]();
  tm0 = new dcomplex[_le]();
  te0 = new dcomplex[_le]();
  vec_at = new dcomplex[_nlemt * _ndm]();
  at = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) at[ai] = vec_at + (ai * _ndm);
}

ParticleDescriptorInclusion::ParticleDescriptorInclusion(const ParticleDescriptorInclusion &rhs) : ParticleDescriptor(rhs) {
  // >>> NEEDED BY CLUSTER AND INCLU <<<
  vec_am0m = new dcomplex[_nlemt * _nlemt];
  np_int nlemts = _nlemt * _nlemt;
  for (np_int ai = 0; ai < nlemts; ai++) vec_am0m[ai] = rhs.vec_am0m[ai];
  vec_fsac = new dcomplex[4];
  vec_sac = new dcomplex[4];
  vec_fsacm = new dcomplex[4];
  for (int vj = 0; vj < 4; vj++) {
    vec_fsac[vj] = rhs.vec_fsac[vj];
    vec_sac[vj] = rhs.vec_sac[vj];
    vec_fsacm[vj] = rhs.vec_fsacm[vj];
  }
  vec_ind3j = new int[(_lm + 1) * _lm];
  np_int vec_ind3j_size = (_lm + 1) * _lm;
  for (np_int ii = 0; ii < vec_ind3j_size; ii++) vec_ind3j[ii] = rhs.vec_ind3j[ii];
  ind3j = new int*[_lm + 1];
  for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
  
  vh = new dcomplex[_ncou * _litpo];
  for (int hi = 0; hi < _ncou * _litpo; hi++) vh[hi] = rhs.vh[hi];
  vj0 = new dcomplex[_nsph * _lmtpo];
  for (int ji = 0; ji < _nsph * _lmtpo; ji++) vj0[ji] = rhs.vj0[ji];
  vyhj = new dcomplex[_ncou * _litpos];
  for (int hj = 0; hj < _ncou * _litpos; hj++) vyhj[hj] = rhs.vyhj[hj];
  vyj0 = new dcomplex[_nsph * _lmtpos];
  for (int jj = 0; jj < _nsph * _lmtpos; jj++) vyj0[jj] = rhs.vyj0[jj];
  am0m = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
  sac = new dcomplex*[2];
  fsac = new dcomplex*[2];
  fsacm = new dcomplex*[2];
  for (int fi = 0; fi < 2; fi++) {
    sac[fi] = vec_sac + (fi * 2);
    fsac[fi] = vec_fsac + (fi * 2);
    fsacm[fi] = vec_fsacm + (fi * 2);
  }
  vintm = new dcomplex[16];
  for (int mi = 0; mi < 16; mi++) vintm[mi] = rhs.vintm[mi];
  scscp = new dcomplex[2];
  ecscp = new dcomplex[2];
  scscpm = new dcomplex[2];
  ecscpm = new dcomplex[2];
  scsc = new double[2];
  ecsc = new double[2];
  scscm = new double[2];
  ecscm = new double[2];
  for (int ci = 0; ci < 2; ci++) {
    scscp[ci] = rhs.scscp[ci];
    ecscp[ci] = rhs.ecscp[ci];
    scscpm[ci] = rhs.scscpm[ci];
    ecscpm[ci] = rhs.ecscpm[ci];
    scsc[ci] = rhs.scsc[ci];
    ecsc[ci] = rhs.ecsc[ci];
    scscm[ci] = rhs.scscm[ci];
    ecscm[ci] = rhs.ecscm[ci];
  }
  v3j0 = new double[_nv3j];
  for (int vj = 0; vj < _nv3j; vj++) v3j0[vj] = rhs.v3j0[vj];
  rac3j = new double[_lmtpo];
  for (int ri = 0; ri < _lmtpo; ri++) rac3j[ri] = rhs.rac3j[ri];
  // >>> NEEDED BY INCLU <<< //
  rm0 = new dcomplex[_le];
  re0 = new dcomplex[_le];
  rmw = new dcomplex[_le];
  rew = new dcomplex[_le];
  tm = new dcomplex[_le];
  te = new dcomplex[_le];
  tm0 = new dcomplex[_le];
  te0 = new dcomplex[_le];
  for (int ti = 0; ti < _le; ti++) {
    rm0[ti] = rhs.rm0[ti];
    re0[ti] = rhs.re0[ti];
    rmw[ti] = rhs.rmw[ti];
    rew[ti] = rhs.rew[ti];
    tm[ti] = rhs.tm[ti];
    te[ti] = rhs.te[ti];
    tm0[ti] = rhs.tm0[ti];
    te0[ti] = rhs.te0[ti];
  }
  vec_at = new dcomplex[_nlemt * _ndm];
  for (int vi = 0; vi < _nlemt * _ndm; vi++) vec_at[vi] = rhs.vec_at[vi];
  at = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) at[ai] = vec_at + (ai * _ndm);
}

#ifdef MPI_VERSION
ParticleDescriptorInclusion::ParticleDescriptorInclusion(const mixMPI *mpidata) : ParticleDescriptor(mpidata) {
  // >>> NEEDED BY CLUSTER AND INCLU <<< //
  vec_am0m = new dcomplex[_nlemt * _nlemt];
  int nlemts = _nlemt * _nlemt;
  MPI_Bcast(vec_am0m, nlemts, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_fsac = new dcomplex[4];
  MPI_Bcast(vec_fsac, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_sac = new dcomplex[4];
  MPI_Bcast(vec_sac, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_fsacm = new dcomplex[4];
  MPI_Bcast(vec_fsacm, 4, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  int vec_ind3j_size = (_lm + 1) * _lm;
  vec_ind3j = new int[vec_ind3j_size];
  MPI_Bcast(vec_ind3j, vec_ind3j_size, MPI_INT, 0, MPI_COMM_WORLD);
  
  vh = new dcomplex[_ncou * _litpo];
  MPI_Bcast(vh, _ncou * _litpo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vj0 = new dcomplex[_nsph * _lmtpo];
  MPI_Bcast(vj0, _nsph * _lmtpo, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vyhj = new dcomplex[_ncou * _litpos];
  MPI_Bcast(vyhj, _ncou * _litpos, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vyj0 = new dcomplex[_nsph * _lmtpos];
  MPI_Bcast(vyj0, _nsph * _lmtpos, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  am0m = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
  sac = new dcomplex*[2];
  fsac = new dcomplex*[2];
  fsacm = new dcomplex*[2];
  for (int fi = 0; fi < 2; fi++) {
    sac[fi] = vec_sac + (fi * 2);
    fsac[fi] = vec_fsac + (fi * 2);
    fsacm[fi] = vec_fsacm + (fi * 2);
  }
  vintm = new dcomplex[16];
  MPI_Bcast(vintm, 16, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  scscp = new dcomplex[2];
  MPI_Bcast(scscp, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  ecscp = new dcomplex[2];
  MPI_Bcast(ecscp, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  scscpm = new dcomplex[2];
  MPI_Bcast(scscpm, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  ecscpm = new dcomplex[2];
  MPI_Bcast(ecscpm, 2, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  scsc = new double[2];
  MPI_Bcast(scsc, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ecsc = new double[2];
  MPI_Bcast(ecsc, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  scscm = new double[2];
  MPI_Bcast(scscm, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ecscm = new double[2];
  MPI_Bcast(ecscm, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  v3j0 = new double[_nv3j];
  MPI_Bcast(v3j0, _nv3j, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ind3j = new int*[_lm + 1];
  for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
  rac3j = new double[_lmtpo];
  MPI_Bcast(rac3j, _lmtpo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // >>> NEEDED BY INCLU <<< //
  rm0 = new dcomplex[_le];
  MPI_Bcast(rm0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  re0 = new dcomplex[_le];
  MPI_Bcast(re0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  rmw = new dcomplex[_le];
  MPI_Bcast(rmw, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  rew = new dcomplex[_le];
  MPI_Bcast(rew, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  tm = new dcomplex[_le];
  MPI_Bcast(tm, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  te = new dcomplex[_le];
  MPI_Bcast(te, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  tm0 = new dcomplex[_le];
  MPI_Bcast(tm0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  te0 = new dcomplex[_le];
  MPI_Bcast(te0, _le, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_at = new dcomplex[_nlemt * _ndm];
  MPI_Bcast(vec_at, _nlemt * _ndm, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  at = new dcomplex*[_nlemt];
  for (int ai = 0; ai < _nlemt; ai++) at[ai] = vec_at + (ai * _ndm);
}
#endif // MPI_VERSION

int ParticleDescriptorInclusion::update_orders(int inner_order, int outer_order) {
  int result = 0;
  bool changed_li = false;
  bool changed_le = false;
  if (inner_order != _li) {
    _li = inner_order;
    changed_li = true;
    _nlim = _li * (_li + 2);
    _litpo = _li + _li + 1;
    _litpos = _litpo * _litpo;
    _ndi = _nsph * _nlim;
    _ndit = 2 * _nsph * _nlim;
    delete[] rmi;
    delete[] vec_rmi;
    delete[] rei;
    delete[] vec_rei;
    vec_rmi = new dcomplex[_li * _nsph]();
    rmi = new dcomplex*[_li];
    vec_rei = new dcomplex[_li * _nsph]();
    rei = new dcomplex*[_li];
    for (int ri = 0; ri < _li; ri++) {
      rmi[ri] = vec_rmi + (_nsph * ri);
      rei[ri] = vec_rei + (_nsph * ri);
    }
    delete[] vh;
    vh = new dcomplex[_ncou * _litpo]();
    delete[] vyhj;
    vyhj = new dcomplex[_ncou * _litpos]();
  }
  if (outer_order != le) {
    _le = outer_order;
    changed_le = true;
    _nlem = _le * (_le + 2);
    _nlemt = 2 * _nlem;
    delete[] w;
    delete[] vec_w;
    vec_w = new dcomplex[_nlemt * 4]();
    w = new dcomplex*[_nlemt];
    for (int wi = 0; wi < _nlemt; wi++) w[wi] = vec_w + (4 * wi);
    delete[] am0m;
    delete[] vec_am0m;
    vec_am0m = new dcomplex[_nlemt * _nlemt]();
    am0m = new dcomplex*[_nlemt];
    for (int ai = 0; ai < _nlemt; ai++) am0m[ai] = vec_am0m + (ai * _nlemt);
    delete[] rm0;
    rm0 = new dcomplex[_le]();
    delete[] re0;
    re0 = new dcomplex[_le]();
    delete[] rmw;
    rmw = new dcomplex[_le]();
    delete[] rew;
    rew = new dcomplex[_le]();
    delete[] tm;
    tm = new dcomplex[_le]();
    delete[] te;
    te = new dcomplex[_le]();
    delete[] tm0;
    tm0 = new dcomplex[_le]();
    delete[] te0;
    te0 = new dcomplex[_le]();
  }
  if (changed_li || changed_le) {
    _lm = (_li > _le) ? _li : _le;
    _lmpo = _lm + 1;
    _lmtpo = _li + _le + 1;
    _lmtpos = _lmtpo * _lmtpo;
    _nv3j = (_lm * (_lm + 1) * (2 * _lm + 7)) / 6;
    _ndm = 2 * (_nsph  * _nlim + _nlem);
    delete[] ind3j;
    delete[] vec_ind3j;
    vec_ind3j = new int[(_lm + 1) * _lm]();
    ind3j = new int*[_lm + 1];
    for (int ii = 0; ii <= _lm; ii++) ind3j[ii] = vec_ind3j + (_lm * ii);
    delete[] vj0;
    vj0 = new dcomplex[_nsph * _lmtpo]();
    delete[] vyj0;
    vyj0 = new dcomplex[_nsph * _lmtpos]();
    delete[] v3j0;
    v3j0 = new double[_nv3j]();
    delete[] rac3j;
    rac3j = new double[_lmtpo]();
    delete[] at;
    delete[] vec_at;
    vec_at = new dcomplex[_nlemt * _ndm]();
    at = new dcomplex*[_nlemt];
    for (int ai = 0; ai < _nlemt; ai++) at[ai] = vec_at + (ai * _ndm);
  }
  return result;
}
// >>> End of ParticleDescriptorInclusion class implementation. <<< //

// >>> ParticleDescriptorSphere class implementation. <<< //
ParticleDescriptorSphere::ParticleDescriptorSphere(GeometryConfiguration *gconf, ScattererConfiguration *sconf) : ParticleDescriptor(gconf, sconf) {
  _class_type = SPHERE_TYPE;
  vec_sas = new dcomplex[4 * _nsph]();
  vec_vints = new dcomplex[16 * _nsph]();

  int nllt = (_nlemt == 0) ? 2 * _nsph * _li * (_li + 2) : _nlemt;
  vec_w = new dcomplex[nllt * 4]();
  w = new dcomplex*[nllt];
  for (int wi = 0; wi < nllt; wi++) w[wi] = vec_w + (4 * wi);
  fsas = new dcomplex[_nsph];
  sscs = new double[_nsph]();
  sexs = new double[_nsph]();
  sabs = new double[_nsph]();
  sqscs = new double[_nsph]();
  sqexs = new double[_nsph]();
  sqabs = new double[_nsph]();
  gcsv = new double[_nsph]();
  sas = new dcomplex**[_nsph];
  vints = new dcomplex*[_nsph];
  for (int vi = 0; vi < _nsph; vi++) {
    vints[vi] = vec_vints + (16 * vi);
    sas[vi] = new dcomplex*[2];
    sas[vi][0] = vec_sas + (4 * vi);
    sas[vi][1] = vec_sas + (4 * vi) + 2;
  }
}

ParticleDescriptorSphere::ParticleDescriptorSphere(const ParticleDescriptorSphere &rhs) : ParticleDescriptor(rhs) {
  // >>> NEEDED BY SPHERE AND CLUSTER <<< //
  vec_sas = new dcomplex[4 * _nsph];
  for (int vsi = 0; vsi < 4 * _nsph; vsi++) vec_sas[vsi] = rhs.vec_sas[vsi];
  vec_vints = new dcomplex[16 * _nsph];
  for (int vvi = 0; vvi < 16 * _nsph; vvi++) vec_vints[vvi] = rhs.vec_vints[vvi];
  sas = new dcomplex**[_nsph];
  vints = new dcomplex*[_nsph];
  for (int vi = 0; vi < _nsph; vi++) {
    vints[vi] = vec_vints + (16 * vi);
    sas[vi] = new dcomplex*[2];
    sas[vi][0] = vec_sas + (4 * vi);
    sas[vi][1] = vec_sas + (4 * vi) + 2;
  }
  
  fsas = new dcomplex[_nsph];
  sscs = new double[_nsph];
  sexs = new double[_nsph];
  sabs = new double[_nsph];
  sqscs = new double[_nsph];
  sqexs = new double[_nsph];
  sqabs = new double[_nsph];
  gcsv = new double[_nsph];
  for (int gi = 0; gi < _nsph; gi++) {
    fsas[gi] = rhs.fsas[gi];
    sscs[gi] = rhs.sscs[gi];
    sexs[gi] = rhs.sexs[gi];
    sabs[gi] = rhs.sabs[gi];
    sqscs[gi] = rhs.sqscs[gi];
    sqexs[gi] = rhs.sqexs[gi];
    sqabs[gi] = rhs.sqabs[gi];
    gcsv[gi] = rhs.gcsv[gi];
  }
}

#ifdef MPI_VERSION
ParticleDescriptorSphere::ParticleDescriptorSphere(const mixMPI *mpidata) : ParticleDescriptor(mpidata) {
  // >>> NEEDED BY SPHERE AND CLUSTER <<< //
  vec_sas = new dcomplex[4 * _nsph];
  MPI_Bcast(vec_sas, 4 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  vec_vints = new dcomplex[16 * _nsph];
  MPI_Bcast(vec_vints, 16 * _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  
  fsas = new dcomplex[_nsph];
  MPI_Bcast(fsas, _nsph, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
  sscs = new double[_nsph];
  MPI_Bcast(sscs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sexs = new double[_nsph];
  MPI_Bcast(sexs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sabs = new double[_nsph];
  MPI_Bcast(sabs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sqscs = new double[_nsph];
  MPI_Bcast(sqscs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sqexs = new double[_nsph];
  MPI_Bcast(sqexs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sqabs = new double[_nsph];
  MPI_Bcast(sqabs, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  gcsv = new double[_nsph];
  MPI_Bcast(gcsv, _nsph, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  sas = new dcomplex**[_nsph];
  vints = new dcomplex*[_nsph];
  for (int vi = 0; vi < nsph; vi++) {
    vints[vi] = vec_vints + (16 * vi);
    sas[vi] = new dcomplex*[2];
    sas[vi][0] = vec_sas + (4 * vi);
    sas[vi][1] = vec_sas + (4 * vi) + 2;
  }
}
#endif // MPI_VERSION

int ParticleDescriptorSphere::update_order(int order) {
  int result = 0;
  if (order != _lm) {
    _lm = order;
    _li = order;
    delete[] rmi;
    delete[] vec_rmi;
    delete[] rei;
    delete[] vec_rei;
    vec_rmi = new dcomplex[_li * _nsph]();
    rmi = new dcomplex*[_li];
    vec_rei = new dcomplex[_li * _nsph]();
    rei = new dcomplex*[_li];
    for (int ri = 0; ri < _li; ri++) {
      rmi[ri] = vec_rmi + (_nsph * ri);
      rei[ri] = vec_rei + (_nsph * ri);
    }
    int nlmmt = 2 * _nsph * _li * (_li + 2);
    delete[] w;
    delete[] vec_w;
    vec_w = new dcomplex[nlmmt * 4]();
    w = new dcomplex*[nlmmt];
    for (int wi = 0; wi < nlmmt; wi++) w[wi] = vec_w + (4 * wi);
  }
  return result;
}
// >>> End of ParticleDescriptorSphere class implementation. <<< //

// >>> ScatteringAngles class implementation. <<< //
ScatteringAngles::ScatteringAngles(GeometryConfiguration *gconf) {
  int isam = gconf->isam;
  _th = gconf->in_theta_start;
  _thstp = gconf->in_theta_step;
  _thlst = gconf->in_theta_end;
  _ths = gconf->sc_theta_start;
  _thsstp = gconf->sc_theta_step;
  _thslst = gconf->sc_theta_end;
  _ph = gconf->in_phi_start;
  _phstp = gconf->in_phi_step;
  _phlst = gconf->in_phi_end;
  _phs = gconf->sc_phi_start;
  _phsstp = gconf->sc_phi_step;
  _phslst = gconf->sc_phi_end;
  double small = 1.0e-3;
  _nth = 0;
  _nph = 0;
  if (_thstp != 0.0) _nth = int((_thlst - _th) / _thstp + small);
  _nth++;
  if (_phstp != 0.0) _nph = int((_phlst - _ph) / _phstp + small);
  _nph++;
  _nths = 0;
  _nphs = 0;
  _thsca = 0.0;
  if (isam > 1) {
    _nths = 1;
    _thsca = _ths - _th;
  } else { // ISAM <= 1
    if (_thsstp == 0.0) _nths = 0;
    else _nths = int ((_thslst - _ths) / _thsstp + small);
    _nths++;
  }
  if (isam >= 1) {
    _nphs = 1;
  } else {
    if (_phsstp == 0.0) _nphs = 0;
    else _nphs = int((_phslst - _phs) / _phsstp + small);
    _nphs++;
  }
  _nk = nth * nph;
  _nks = nths * nphs;
  _nkks = nk * nks;
}

ScatteringAngles::ScatteringAngles(const ScatteringAngles &rhs) {
  _th = rhs._th;
  _thstp = rhs._thstp;
  _thlst = rhs._thlst;
  _ths = rhs._ths;
  _thsstp = rhs._thsstp;
  _thslst = rhs._thslst;
  _ph = rhs._ph;
  _phstp = rhs._phstp;
  _phlst = rhs._phlst;
  _phs = rhs._phs;
  _phsstp = rhs._phsstp;
  _phslst = rhs._phslst;
  _nth = rhs._nth;
  _nph = rhs._nph;
  _nths = rhs._nths;
  _nphs = rhs._nphs;
  _thsca = rhs._thsca;
  _nk = rhs._nk;
  _nks = rhs._nks;
  _nkks = rhs._nkks;
}

#ifdef MPI_VERSION
ScatteringAngles::ScatteringAngles(const mixMPI *mpidata) {
  MPI_Bcast(&_th, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thlst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ths, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thsstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thslst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ph, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phlst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phsstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phslst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nth, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nph, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nths, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nphs, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thsca, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nks, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nkks, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void ScatteringAngles::mpibcast(const mixMPI *mpidata) {
  MPI_Bcast(&_th, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thlst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ths, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thsstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thslst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_ph, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phlst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phsstp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_phslst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nth, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nph, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nths, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nphs, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_thsca, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nks, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&_nkks, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
#endif
// >>> End of ScatteringAngles class implementation. <<< //
