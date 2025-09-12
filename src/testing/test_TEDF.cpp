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

//! \file test_TEDF.cpp

#include <cstdio>
#include <exception>
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

#ifndef INCLUDE_FILE_IO_H_
#include "../include/file_io.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

using namespace std;

/*! \brief Main program execution body.
 *
 * This program executes a test to compare whether three configuration
 * instances, loaded respectively from an EDFB configuration file, a
 * legacy binary, or a HDF5 binary are actually equivalent. The test
 * writes a result message to `stdout` then it returns 0 (OS flag for
 * a successful process) or some kind of error code, depending on
 * whether the test files were found all equal or not. The test accepts
 * three command line arguments: the name of the EDFB configuration
 * file, the name of the legacy binary file and the name of the HDF5
 * binary configuration file.
 *
 * \param argc: `int` Number of command line arguments
 * \param argv: `char **` Array of command argument character strings.
 * \return result: `int` Can be: 0 (all files equal); 1 (EDFB and
 * legacy binary are different); 10 (EDFB and HDF5 are different);
 * 100 (legacy and HDF5 are different). In case more differences are
 * found, the error codes sum up together (e.g. 111 means all files
 * are different).
 */
int main(int argc, char **argv) {
  int result = 0;
  string dedfb_file = "DEDFB";
  string legacy_file = "c_TEDF";
  string hdf5_file = "c_TEDF.hd5";
  if (argc == 4) {
    dedfb_file = string(argv[1]);
    legacy_file = string(argv[2]);
    hdf5_file = string(argv[3]);
  }
  ScattererConfiguration *a = NULL, *b = NULL, *c = NULL;
  try { a = ScattererConfiguration::from_dedfb(dedfb_file); }
  catch (...) {
    printf("Failed to open legacy configuration file.\n");
    return 2;
  }
  try { b = ScattererConfiguration::from_binary(legacy_file); }
  catch (...) {
    printf("Failed to open legacy binary file.\n");
    return 20;
  }
  try { c = ScattererConfiguration::from_binary(hdf5_file, "HDF5"); }
  catch (...) {
    printf("Failed to open HDF5 configuration file.\n");
    return 200;
  }
  if (*a == *b) printf("Configuration objects a and b are equal.\n");
  else {
    printf("Configuration objects a and b are different.\n");
    result += 1;
  }
  if (*a == *c) printf("Configuration objects a and c are equal.\n");
  else {
    printf("Configuration objects a and c are different.\n");
    result += 10;
  }
  if (*c == *b) printf("Configuration objects c and b are equal.\n");
  else {
    printf("Configuration objects c and b are different.\n");
    result += 100;
  }
  if (a != NULL) delete a;
  if (b != NULL) delete b;
  if (c != NULL) delete c;
  return result;
}
