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

/*! \file np_trapping.cpp
 *
 * \brief Trapping problem handler.
 */
#include <chrono>
#include <cstdio>
#include <string>

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

using namespace std;

extern void frfme(string data_file, string output_path);
extern void lffft(string data_file, string output_path);

/*! \brief Main program execution body.
 *
 * \param argc: `int`
 * \param argv: `char **`
 * \return result: `int`
 */
int main(int argc, char **argv) {
  chrono::time_point<chrono::high_resolution_clock> t_start = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed;
  string frfme_data_file = "../../test_data/trapping/DFRFME";
  string lffft_data_file = "../../test_data/trapping/DLFFFT";
  string output_path = ".";
  string message;
  Logger logger(LOG_DEBG);
  if (argc == 4) {
    frfme_data_file = string(argv[1]);
    lffft_data_file = string(argv[2]);
    output_path = string(argv[3]);
  }
  frfme(frfme_data_file, output_path);
  lffft(lffft_data_file, output_path);
  elapsed = chrono::high_resolution_clock::now() - t_start;
  message = "INFO: calculation lasted " + to_string(elapsed.count()) + "s.\n";
  logger.log(message);
  return 0;
}
