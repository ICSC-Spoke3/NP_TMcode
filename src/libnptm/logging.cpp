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
 * \ingroup support
 */
#include <cstdio>
#include <string>

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

using namespace std;

//! \cond
// Tags of the logging levels, indexed by LOG_DEBG, LOG_INFO, LOG_WARN, LOG_ERRO.
static const char *level_tags[4] = {"[DEBUG] ", "[INFO] ", "[WARN] ", "[ERROR] "};

static const char *tag_of(int level) {
  if (level < LOG_DEBG) level = LOG_DEBG;
  if (level > LOG_ERRO) level = LOG_ERRO;
  return level_tags[level];
}

// Whether text starts with a level tag at position pos.
static bool starts_with_tag(const string& text, size_t pos) {
  for (int ti = 0; ti < 4; ti++) {
    if (text.compare(pos, string(level_tags[ti]).size(), level_tags[ti]) == 0) return true;
  }
  return false;
}
//! \endcond

Logger::Logger(int threshold, FILE *logging_output, FILE *error_output) {
  last_message = new string("");
  log_threshold = threshold;
  log_output = logging_output;
  err_output = error_output;
  repetitions = 0;
  log_at_line_start = true;
  err_at_line_start = true;
}

Logger::~Logger() {
  delete last_message;
}

string Logger::tag_lines(const string& message, const char *tag, bool &at_line_start) const {
  string result;
  result.reserve(message.size() + 16);
  for (size_t ci = 0; ci < message.size(); ci++) {
    const char c = message[ci];
    if (at_line_start && c != '\n') {
      if (!starts_with_tag(message, ci)) result += tag;
      at_line_start = false;
    }
    result += c;
    if (c == '\n') at_line_start = true;
  }
  return result;
}

void Logger::err(const std::string& message) const {
  // One thread at a time, so that the lines of different messages are not
  // mixed and the line-start state stays consistent.
#pragma omp critical(np_logger_output)
  {
    string text = tag_lines(message, tag_of(LOG_ERRO), err_at_line_start);
    fprintf(err_output, "%s", text.c_str());
    fflush(err_output);
  }
}

void Logger::flush(int level) {
  string summary = "\"" + *last_message + "\" issued " + to_string(repetitions);
  if (repetitions == 1) summary += " time.\n";
  else summary += " times.\n";
  log(summary, level);
  delete last_message;
  last_message = new string("");
  repetitions = 0;
}

void Logger::log(const std::string& message, int level) const {
  if (level == LOG_ERRO) err(message);
  else {
    if (level >= log_threshold) {
#pragma omp critical(np_logger_output)
      {
	string text = tag_lines(message, tag_of(level), log_at_line_start);
	fprintf(log_output, "%s", text.c_str());
	fflush(log_output);
      }
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
