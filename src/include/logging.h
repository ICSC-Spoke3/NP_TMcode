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

/*! \file logging.h
 *
 * \brief Definition of the logging system.
 * \ingroup support
 */
#ifndef INCLUDE_LOGGING_H_
#define INCLUDE_LOGGING_H_

/*! \addtogroup support
 * @{
 */

//! \cond
#define NPTM_STRINGIFY_(x) #x
#define NPTM_STRINGIFY(x) NPTM_STRINGIFY_(x)
//! \endcond

#ifdef NPTM_VERSION_ID
/*! \brief Code version, from the "version" key of `CITATION.cff`.
 *
 * `configure` reads the key and defines `NPTM_VERSION_ID` (without
 * quotes); a build without it reports "unknown". Rerun `configure` after
 * changing the version.
 */
#define NPTM_VERSION NPTM_STRINGIFY(NPTM_VERSION_ID)
#else
//! \brief Code version: unknown, because `configure` did not define `NPTM_VERSION_ID`.
#define NPTM_VERSION "unknown"
#endif

//! \brief Debug level logging (maximum verbosity).
#define LOG_DEBG 0
//! \brief Standard information level logging (default).
#define LOG_INFO 1
//! \brief Warning level logging (almost quiet).
#define LOG_WARN 2
//! \brief Error level logging (silent, unless in pain).
#define LOG_ERRO 3

//! \brief Macro to stringize code lines
#define TOSTRING(ARG) string(#ARG)

/*! \brief Logger class.
 *
 * Loggers are objects used to track the execution of a code, reporting activities
 * such as function calls and parameter settings. They can be used to inform the
 * user about the execution of operations at runtime (e.g. by printing messages to
 * the terminal), as well as to record the execution history in appropriate log
 * files. The `Logger` class offers an implementation of logging system complying
 * with the requirements of the NP_TMcode project.
 *
 * The Logger class is designed to work with open files. It is a user responsibility
 * to check that the required log files are properly opened before use, and closed
 * thereafter, if they are not the standard `stdout` and `stderr` streams.
 *
 * Every line written by a logger starts with a tag that names its level:
 * `[DEBUG] `, `[INFO] `, `[WARN] ` or `[ERROR] `. The logger adds the tag, so
 * messages must not include it. A message may be written in several calls:
 * the tag is added only at the start of a line (a call that continues a line
 * gets none), and a line that already starts with a tag (e.g. the text of an
 * exception carrying `[ERROR] `) is not tagged again.
 */
class Logger {
 protected:
  //! \brief Pointer to error stream.
  FILE *err_output;
  //! \brief Pointer to logging stream.
  FILE *log_output;
  //! \brief Last logged message.
  std::string *last_message;
  //! \brief Threshold of logging level.
  int log_threshold;
  //! \brief Number of identical message repetitions.
  long repetitions;
  //! \brief Whether the next character of the logging stream starts a line.
  mutable bool log_at_line_start;
  //! \brief Whether the next character of the error stream starts a line.
  mutable bool err_at_line_start;

  /*! \brief Add a tag at the start of each line of a message.
   *
   * \param message: `const string &` The message to be tagged.
   * \param tag: `const char *` The tag to add.
   * \param at_line_start: `bool &` Whether the message starts a line
   * (updated for the next message).
   * \return tagged: `string` The tagged message.
   */
  std::string tag_lines(const std::string& message, const char *tag, bool &at_line_start) const;

 public:
  /*! \brief Logger instance constructor.
   *
   * \param threshold: `int` Threshold of the messages to be included in log. Can
   * be `LOG_DEBG` (log everything), `LOG_INFO` (give detailed information),
   * `LOG_WARN` (log odd looking effects), or `LOG_ERRO` (print error messages,
   * always active). The default behaviour is `LOG_WARN`.
   * \param logging_output: `FILE *` Pointer to an open output file for common messages
   * (optional, default is `stdout`).
   * \param error_output: `FILE *` Pointer to an open output file for error messages
   * (optional, default is `stderr`).
   */
  Logger(int threshold, FILE *logging_output = stdout, FILE *error_output = stderr);

  /*! \brief Logger instance destroyer.
   */
  ~Logger();

  /*! \brief Print a message to the error output, tagged as `[ERROR] `.
   *
   * \param message: `string` The message to be printed.
   */
  void err(const std::string& message) const;

  /*! \brief Print a summary of recurrent messages and clear the stack.
   *
   * \param level: `int` The priority level (default is `LOG_DEBG = 0`).
   */
  void flush(int level=LOG_DEBG);

  /*! \brief Print a message, depending on its logging level.
   *
   * The message is tagged with its level. `LOG_ERRO` messages go to the
   * error output, like those printed with `err()`.
   *
   * \param message: `string` The message to be printed.
   * \param level: `int` The priority level (default is `LOG_INFO = 1`).
   */
  void log(const std::string& message, int level=LOG_INFO) const;
  
  /*! \brief Push a recurrent message to the message stack.
   *
   * When a long stream of identical messages is expected, it may be more
   * convenient to put them in a stack. The stack is the error output stream,
   * so that no memory is consumed. After the call stack is over, or is the
   * message changes, the stack is flushed, meaning that a summary message
   * is written to the logging output and the stack counter is reset to 0.
   *
   * \param message: `string` The message to be stacked.
   */
  void push(const std::string& message);
};

//! @}
#endif
