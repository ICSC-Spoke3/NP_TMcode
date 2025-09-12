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

//! \file test_TTMS.cpp

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

#ifndef INCLUDE_TRANSITIONMATRIX_H_
#include "../include/TransitionMatrix.h"
#endif

using namespace std;

/*! \brief Main program execution body.
 *
 * This program executes a test to compare whether two transition
 * matrix instances, loaded respectively from a legacy and a HDF5
 * binary file are actually equivalent. The test writes a result
 * message to `stdout` then it returns 0 (OS flag for a successful
 * process) or 1 (OS flag for failing process) depending on whether
 * the two instances were considered equivalent or not.
 *
 * \param argc: `int` Number of command line arguments
 * \param argv: `char **` Array of command argument character strings.
 * \return result: `int` Can be 0 (files are equal) or 1 (files are
 * different).
 */
int main(int argc, char **argv) {
  int result = 0;
  TransitionMatrix *a, *b;
  string legacy_file = "c_TTMS";
  string hdf5_file = "c_TTMS.hd5";
  if (argc == 3) {
    legacy_file = string(argv[1]);
    hdf5_file = string(argv[2]);
  }
  a = TransitionMatrix::from_binary(legacy_file);
  b = TransitionMatrix::from_binary(hdf5_file, "HDF5");
  if (*a == *b) printf("Transition matrices a and b are equal.\n");
  else {
    printf("Transition matrices a and b are different.\n");
    result = 1;
  }
  return result;
}
