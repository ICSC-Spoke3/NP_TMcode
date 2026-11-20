/*! \file Parsers.h
 *
 * \brief A library of functions designed to parse formatted input
 * into memory.
 * \ingroup io
 */

#ifndef INCLUDE_PARSERS_H_
#define INCLUDE_PARSERS_H_

#include <string>

#ifndef INCLUDE_ERRORS_H_
#include "errors.h"
#endif

/*! \addtogroup io
 * @{
 */

/*! \brief Load a text file as a sequence of strings in memory.
 *
 * The configuration of the field expansion code in FORTRAN uses
 * shared memory access and file I/O operations managed by different
 * functions. Although this approach could be theoretically replicated,
 * it is more convenient to handle input and output to distinct files
 * using specific functions. load_file() helps in the task of handling
 * input such as configuration files or text data structures that need
 * to be loaded entirely. The function performs a line-by line scan of
 * the input file and returns an array of strings that can be later
 * parsed and ingested by the concerned code blocks. An optional pointer
 * to integer allows the function to keep track of the number of file
 * lines that were read, if needed.
 *
 * \param file_name: `string` The path of the file to be read.
 * \param count: `int*` Pointer to an integer recording the number of
 * read lines [OPTIONAL, default=NULL].
 * \return array_lines `string*` An array of strings, one for each input
 * file line.
 * \throws InputError If the file cannot be opened.
 */
std::string *load_file(const std::string& file_name, int *count);

/*! \brief Run a parsing function, reporting malformed data as a `ParseError`.
 *
 * The parsers of the input files convert text with `std::stoi()`,
 * `std::stod()` and similar functions, which throw `std::invalid_argument`
 * or `std::out_of_range` on malformed text. This function runs a parsing
 * body and turns those exceptions into a `ParseError` that names the file
 * and the line being read, so that the user can find the problem.
 *
 * \param file_name: `const string &` Name of the file being parsed.
 * \param current_line: `const int &` Index (from 0) of the line being
 * parsed, kept up to date by the parsing body.
 * \param body: Function object with the parsing code, called with no
 * arguments.
 * \return result: The value returned by `body`.
 * \throws ParseError If `body` throws `std::invalid_argument` or `std::out_of_range`.
 */
template <class F> auto parse_guarded(const std::string& file_name, const int& current_line, F body) -> decltype(body()) {
  try {
    return body();
  } catch (const std::invalid_argument &ex) {
    throw ParseError(
      "malformed value in \"" + file_name + "\" at line " + std::to_string(current_line + 1)
    );
  } catch (const std::out_of_range &ex) {
    throw ParseError(
      "value out of range in \"" + file_name + "\" at line " + std::to_string(current_line + 1)
    );
  }
}

/*! \brief Access a line of a loaded file, checking that it exists.
 *
 * \param file_lines: `string *` The lines returned by `load_file()`.
 * \param num_lines: `int` The number of lines.
 * \param index: `int` Index (from 0) of the requested line.
 * \param file_name: `const string &` Name of the file, for the error message.
 * \param current_line: `int &` Updated to `index`, for later error messages.
 * \return line: `string &` The requested line.
 * \throws ParseError If the file has no such line (it ends too early).
 */
inline std::string& checked_line(
  std::string *file_lines, int num_lines, int index, const std::string& file_name, int &current_line
) {
  if (index < 0 || index >= num_lines) {
    throw ParseError(
      "unexpected end of file \"" + file_name + "\" (" + std::to_string(num_lines)
      + " lines, line " + std::to_string(index + 1) + " needed)"
    );
  }
  current_line = index;
  return file_lines[index];
}

//! @}
#endif /* INCLUDE_PARSERS_H_ */
