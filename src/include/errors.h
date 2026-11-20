/*! \file errors.h
 *
 * \brief Collection of proprietary code exceptions.
 *
 * There are many circumstances that can prevent the correct execution
 * of a code. These range from user mistakes, to improper configuration,
 * to unsupported hardware and all the way up to various system failures.
 * Although it is not possible to grant proper execution in all cases,
 * it is often possible to design a code in such a way that the program
 * detects unexpected conditions, informs the user and takes the proper
 * actions, eventually stopping without crash, if no other options are
 * available. C++ handles such unexpected circumstances by means of
 * `exceptions`. These are special procedures that can be launched
 * whenever an unexpected situation occurs and they allow to restore the
 * code work-flow and attempt recovery. Exceptions can be divided in
 * different cathegories, which respond to various types of problems.
 * This library contains a set of exceptions designed to the most common
 * problems that may occur while executing an application of the `NP_TMcode`
 * suite.
 *
 * All the exceptions of the suite derive from `NptmError`, which carries a
 * description of the problem and the exit status that a program returns when
 * the exception stops it. The exit status identifies the category of the
 * problem (see `NptmExitCode`), so that scripts can react to it:
 *
 * | Category | Class | Exit status |
 * | --- | --- | --- |
 * | Internal error | `NptmError` | 1 |
 * | Usage error | `UsageError` | 2 |
 * | Input file missing or unreadable | `InputError` | 3 |
 * | Input data malformed or inconsistent | `ParseError` | 4 |
 * | Output folder or file not writable | `OutputError` | 5 |
 * | Not enough memory | `MemoryError` | 6 |
 * | GPU or accelerated library failure | `DeviceError` | 7 |
 *
 * Library code reports problems by throwing these exceptions and never
 * stops the program itself: each program catches them in `main()`, logs the
 * description and returns the exit status.
 * \ingroup support
 */

#ifndef INCLUDE_ERRORS_H_
#define INCLUDE_ERRORS_H_

#include <stdexcept>
#include <string>

/*! \addtogroup support
 * @{
 */

/*! \brief Exit status of the programs, by category of problem.
 */
enum NptmExitCode {
  NPTM_EXIT_SUCCESS = 0,  //!< Successful execution.
  NPTM_EXIT_INTERNAL = 1, //!< Internal error (a bug, or an unexpected condition).
  NPTM_EXIT_USAGE = 2,    //!< Wrong command line.
  NPTM_EXIT_INPUT = 3,    //!< Input file missing or unreadable.
  NPTM_EXIT_PARSE = 4,    //!< Input data malformed or inconsistent.
  NPTM_EXIT_OUTPUT = 5,   //!< Output folder or file not writable.
  NPTM_EXIT_MEMORY = 6,   //!< Not enough host or device memory.
  NPTM_EXIT_DEVICE = 7    //!< GPU initialization or accelerated library failure.
};

/*! \brief Base class of the exceptions of the suite.
 *
 * It describes the problem and carries the exit status of its category. Used
 * directly, it reports an internal error.
 */
class NptmError: public std::runtime_error {
protected:
  //! \brief Exit status of the category of the problem.
  int _exit_code;
  //! \brief Whether all the MPI processes know about the problem.
  bool _all_processes;

public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   * \param exit_code: `int` Exit status (optional, default is NPTM_EXIT_INTERNAL).
   */
  NptmError(const std::string& message, int exit_code = NPTM_EXIT_INTERNAL) :
    std::runtime_error(message), _exit_code(exit_code), _all_processes(false) { }

  /**
   * \brief Exit status of the category of the problem.
   *
   * \return code: `int` The exit status.
   */
  int exit_code() const noexcept { return _exit_code; }

  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  virtual const char* category() const noexcept { return "internal error"; }

  /**
   * \brief Whether all the MPI processes know about the problem.
   *
   * Problems detected at a point where the MPI processes agree on their
   * status stop all the processes in the same way, so they can finalize MPI
   * normally. The others require `MPI_Abort()`.
   *
   * \return all_processes: `bool` True if all the processes know.
   */
  bool all_processes() const noexcept { return _all_processes; }

  /**
   * \brief Mark the problem as known to all the MPI processes.
   */
  void set_all_processes() noexcept { _all_processes = true; }
};

/*! \brief Exception for wrong command lines.
 */
class UsageError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   */
  UsageError(const std::string& message) : NptmError(message, NPTM_EXIT_USAGE) { }
  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  const char* category() const noexcept override { return "usage error"; }
};

/*! \brief Exception for missing or unreadable input files.
 */
class InputError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   */
  InputError(const std::string& message) : NptmError(message, NPTM_EXIT_INPUT) { }
  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  const char* category() const noexcept override { return "input error"; }
};

/*! \brief Exception for malformed or inconsistent input data.
 */
class ParseError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   */
  ParseError(const std::string& message) : NptmError(message, NPTM_EXIT_PARSE) { }
  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  const char* category() const noexcept override { return "input data error"; }
};

/*! \brief Exception for output folders or files that cannot be written.
 */
class OutputError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   */
  OutputError(const std::string& message) : NptmError(message, NPTM_EXIT_OUTPUT) { }
  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  const char* category() const noexcept override { return "output error"; }
};

/*! \brief Exception for calculations that do not fit in the available memory.
 */
class MemoryError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   */
  MemoryError(const std::string& message) : NptmError(message, NPTM_EXIT_MEMORY) { }
  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  const char* category() const noexcept override { return "memory error"; }
};

/*! \brief Exception for GPU and accelerated library failures.
 */
class DeviceError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param message: `const string &` Description of the problem.
   */
  DeviceError(const std::string& message) : NptmError(message, NPTM_EXIT_DEVICE) { }
  /**
   * \brief Name of the category of the problem.
   *
   * \return category: `const char *` The category name.
   */
  const char* category() const noexcept override { return "device error"; }
};

/*! \brief Exception for a process that stops because another one failed.
 *
 * When an MPI process detects a problem at a point where the processes agree
 * on their status, the others stop with this exception, which carries the
 * exit status of the problem. Its message is not shown to the user, because
 * the process that detected the problem reports it.
 */
class PeerProcessError: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param exit_code: `int` Exit status of the problem in the other process.
   */
  PeerProcessError(int exit_code) :
    NptmError("another MPI process stopped with exit status " + std::to_string(exit_code), exit_code) {
    _all_processes = true;
  }
};

/*! \brief Exception for wrong OutputInfo NULL calls.
 */
class UnrecognizedOutputInfo: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param requested: `int` The index that was requested.
   */
  UnrecognizedOutputInfo(int requested) :
    NptmError("passed parameter " + std::to_string(requested) + ", but only 1 is allowed") { }
};

/*! \brief Exception for out of bounds List requests.
 */
class ListOutOfBoundsException: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param requested: `int` The index that was requested.
   * \param min: `int` The minimum index allowed by the list.
   * \param max: `int` The maximum index allowed by the list.
   */
  ListOutOfBoundsException(int requested, int min, int max) :
    NptmError(
      "requested index " + std::to_string(requested) + " out of list allowed bounds ["
      + std::to_string(min) + ", " + std::to_string(max - 1) + "]"
    ) { }
};

/*! \brief Exception for data structures that cannot be allocated.
 */
class ObjectAllocationException: public MemoryError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `const string &` Description of the problem.
   */
  ObjectAllocationException(const std::string& problem) : MemoryError(problem) { }
};

/*! \brief Exception for configuration files that cannot be opened.
 */
class OpenConfigurationFileException: public InputError {
protected:
  //! \brief Name of the file that was accessed.
  std::string _file_name;

public:
  /**
   * \brief Exception instance constructor.
   *
   * \param name: `const string &` Name of the file that was accessed.
   */
  OpenConfigurationFileException(const std::string& name) :
    InputError("cannot open configuration file \"" + name + "\""), _file_name(name) { }

  /**
   * \brief Name of the file that could not be opened.
   *
   * \return name: `const string &` The file name.
   */
  const std::string& file_name() const noexcept { return _file_name; }
};

/*! \brief Exception for access requests out of matrix bounds.
 */
class MatrixOutOfBoundsException: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `const string &` Description of the problem that occurred.
   */
  MatrixOutOfBoundsException(const std::string& problem) : NptmError(problem) { }
};

/*! \brief Exception for unrecognized configuration data sets.
 */
class UnrecognizedConfigurationException: public ParseError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `const string &` Description of the problem that occurred.
   */
  UnrecognizedConfigurationException(const std::string& problem) : ParseError(problem) { }
};

/*! \brief Exception for unrecognized file formats.
 */
class UnrecognizedFormatException: public ParseError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `const string &` Description of the problem that occurred.
   */
  UnrecognizedFormatException(const std::string& problem) : ParseError(problem) { }
};

/*! \brief Exception for unrecognized parameters.
 */
class UnrecognizedParameterException: public NptmError {
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `const string &` Description of the problem that occurred.
   */
  UnrecognizedParameterException(const std::string& problem) : NptmError(problem) { }
};

//! @}
#endif
