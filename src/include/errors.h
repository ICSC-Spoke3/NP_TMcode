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
 */

#ifndef INCLUDE_ERRORS_H_
#define INCLUDE_ERRORS_H_

/*! \brief Exception for wrong OutputInfo NULL calls.
 */
class UnrecognizedOutputInfo: public std::exception {
protected:
  //! Description of the problem.
  std::string message;

public:
  /**
   * \brief Exception instance constructor.
   *
   * \param requested: `int` The index that was requested.
   */
  UnrecognizedOutputInfo(int requested) {
    message = "Error: passed parameter " + std::to_string(requested)
      + ", but only 1 is allowed";
  }
  
  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return message.c_str();
  }
};

/*! \brief Exception for out of bounds List requests.
 */
class ListOutOfBoundsException: public std::exception {
protected:
  //! Description of the problem.
  std::string message;

public:
  /**
   * \brief Exception instance constructor.
   *
   * \param requested: `int` The index that was requested.
   * \param min: `int` The minimum index allowed by the list.
   * \param max: `int` The maximum index allowed by the list.
   */
  ListOutOfBoundsException(int requested, int min, int max) {
    message = "Error: requested index " + std::to_string(requested)
      + " out of list allowed bounds [" + std::to_string(min) + ", "
      + std::to_string(max - 1) + "]";
  }
  
  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return message.c_str();
  }
};

/*! \brief Exception for open file error handlers.
 */
class OpenConfigurationFileException: public std::exception {
protected:
  //! \brief Name of the file that was accessed.
  std::string file_name;

public:
  /**
   * \brief Exception instance constructor.
   *
   * \param name: `string` Name of the file that was accessed.
   */
  OpenConfigurationFileException(const std::string& name) { file_name = name; }

  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return file_name.c_str();
  }
};

/*! \brief Exception for access requests out of matrix bounds.
 */
class MatrixOutOfBoundsException: public std::exception {
protected:
  //! Description of the problem.
  std::string message;
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `string` Description of the problem that occurred.
   */
  MatrixOutOfBoundsException(const std::string& problem) { message = problem; }
  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return message.c_str();
  }
};

/*! \brief Exception for unrecognized configuration data sets.
 */
class UnrecognizedConfigurationException: public std::exception {
protected:
  //! Description of the problem.
  std::string message;
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `string` Description of the problem that occurred.
   */
  UnrecognizedConfigurationException(const std::string& problem) { message = problem; }
  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return message.c_str();
  }
};

/*! \brief Exception for unrecognized file formats.
 */
class UnrecognizedFormatException: public std::exception {
protected:
  //! Description of the problem.
  std::string message;
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `string` Description of the problem that occurred.
   */
  UnrecognizedFormatException(const std::string& problem) { message = problem; }
  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return message.c_str();
  }
};

/*! \brief Exception for unrecognized parameters.
 */
class UnrecognizedParameterException: public std::exception {
protected:
  //! Description of the problem.
  std::string message;
public:
  /**
   * \brief Exception instance constructor.
   *
   * \param problem: `string` Description of the problem that occurred.
   */
  UnrecognizedParameterException(const std::string& problem) { message = problem; }
  /**
   * \brief Exception message.
   */
  virtual const char* what() const throw() {
    return message.c_str();
  }
};

#endif
