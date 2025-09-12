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

/*! \file np_cluster.cpp
 *
 * \brief Cluster of spheres scattering problem handler.
 *
 * This program emulates the execution work-flow originally handled by the
 * FORTRAN EDFB and CLU codes, which undertook the task of solving the
 * scattering calculation for the case of a cluster of spheres, with one or
 * more materials. The program is designed to work in two modes: emulation and
 * enhanced. The emulation mode is activated by invoking the program without
 * arguments. In this case, the code looks for hard-coded default input and
 * writes output in the execution folder, replicating the behaviour of the
 * original FORTRAN code. Advanced mode, instead, is activated by passing
 * command line arguments that locate the desired input files and a valid
 * folder to write the output into. The emulation mode is useful for testing,
 * while the advanced mode implements the possibility to change input and
 * output options, without having to modify the code.
 */

#include <cstdio>
#include <string>

#ifdef USE_MPI
#ifndef MPI_VERSION
#include <mpi.h>
#endif
#endif

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

using namespace std;

extern void cluster(const string& config_file, const string& data_file, const string& output_path, const mixMPI *mpidata);

/*! \brief Main program entry point.
 *
 * This is the starting point of the execution flow. Here we may choose
 * how to configure the code, e.g. by loading a legacy configuration file
 * or some otherwise formatted configuration data set. The code can be
 * linked to a luncher script or to a GUI oriented application that performs
 * the configuration and runs the main program.
 *
 * \param argc: `int` The number of arguments given in command-line.
 * \param argv: `char **` The vector of command-line arguments.
 * \return result: `int` An exit code passed to the OS (0 for succesful execution).
 */
int main(int argc, char **argv) {
  int ierr = 0;
#ifdef MPI_VERSION
  ierr = MPI_Init(&argc, &argv);
  // create and initialise class with essential MPI data
  mixMPI *mpidata = new mixMPI(MPI_COMM_WORLD);
#else
  // create a the class with dummy data if we are not using MPI at all
  mixMPI *mpidata = new mixMPI();
#endif
  string config_file = "../../test_data/cluster/DEDFB";
  string data_file = "../../test_data/cluster/DCLU";
  string output_path = ".";
  if (argc == 4) {
    config_file = string(argv[1]);
    data_file = string(argv[2]);
    output_path = string(argv[3]);
  }
  cluster(config_file, data_file, output_path, mpidata);
#ifdef MPI_VERSION
  MPI_Finalize();
#endif
  delete mpidata;
  return ierr;
}
