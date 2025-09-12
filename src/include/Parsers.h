/*! \file Parsers.h
 *
 * \brief A library of functions designed to parse formatted input
 * into memory.
 */

#ifndef FILE_NOT_FOUND_ERROR
//! Error code if a file is not found.
#define FILE_NOT_FOUND_ERROR 21
#endif

#ifndef INCLUDE_PARSERS_H_
#define INCLUDE_PARSERS_H_

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
 */
std::string *load_file(const std::string& file_name, int *count);

#endif /* INCLUDE_PARSERS_H_ */
