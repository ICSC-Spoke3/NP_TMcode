/*! \file Parsers.cpp
 *
 * \brief Implementation of the parsing functions.
 */

#include <exception>
#include <fstream>
#include <string>

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_PARSERS_H_
#include "../include/Parsers.h"
#endif

std::string *load_file(const std::string& file_name, int *count = 0) {
  std::fstream input_file(file_name.c_str(), std::ios::in);
  List<std::string> *file_lines = new List<std::string>();
  std::string line;
  if (input_file.is_open()) {
    getline(input_file, line);
    file_lines->set(0, line);
    while (getline(input_file, line)) {
      file_lines->append(line);
    }
    input_file.close();
  } else {
	  throw FILE_NOT_FOUND_ERROR;
  }
  std::string *array_lines = file_lines->to_array();
  if (count != 0) *count = file_lines->length();
  delete file_lines;
  return array_lines;
}
