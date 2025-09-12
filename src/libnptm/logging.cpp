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

/*! \file logging.cpp
 *
 * \brief Implementation of the logging system.
 */
#include <cstdio>
#include <string>

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

using namespace std;

Logger::Logger(int threshold, FILE *logging_output, FILE *error_output) {
  last_message = new string("");
  log_threshold = threshold;
  log_output = logging_output;
  err_output = error_output;
  repetitions = 0;
}

Logger::~Logger() {
  delete last_message;
}

void Logger::err(const std::string& message) {
  fprintf(err_output, "%s", message.c_str());
  fflush(err_output);
}

void Logger::flush(int level) {
  string summary = "\"" + *last_message + "\" issued " + to_string(repetitions);
  if (repetitions == 1) summary += " time.\n";
  else summary += " times.\n";
  if (level == LOG_ERRO) err(summary);
  else {
    if (level >= log_threshold) {
      fprintf(log_output, "%s", summary.c_str());
      fflush(log_output);
    }
  }
  delete last_message;
  last_message = new string("");
  repetitions = 0;
}

void Logger::log(const std::string& message, int level) {
  if (level == LOG_ERRO) err(message);
  else {
    if (level >= log_threshold) {
      fprintf(log_output, "%s", message.c_str());
      fflush(log_output);
    }
  }
}

void Logger::push(const std::string& message) {
  if (repetitions > 0) {
    if (message.compare(*last_message) != 0) {
      flush(LOG_DEBG);
    }
  }
  log(message, LOG_DEBG);
  delete last_message;
  last_message = new string(message);
  repetitions++;
}
