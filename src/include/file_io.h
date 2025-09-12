/*! \file file_io.h
 *
 * \brief Library to handle I/O operations with files.
 */
#ifndef INCLUDE_FILE_IO_H_
#define INCLUDE_FILE_IO_H_

#include <vector>

class mixMPI;

/*! \class FileSchema
 *
 * \brief File content descriptor.
 *
 * Accessing binary files requires detailed knowledge of their contents. The `FileSchema`
 * class is intended to encapsulate this information and use it as a wrapper to control
 * I/O operations towards different file formats. Any file can be thought of as a sequence
 * of records, which may further contain arbitrarily complex structures. By describing the
 * structure of records, it is possible to support virtually any type of format.
 */
class FileSchema {
 protected:
  //! \brief Number of records conained in the file.
  int num_records;
  //! \brief Array of record names.
  std::string *record_names;
  //! \brief Array of record descriptors.
  std::string *record_types;

 public:
  /*! \brief FileSchema instance constructor.
   *
   * \param num_rec: `int` Number of records in the file.
   * \param rec_types: `string *` Description of the records in the file.
   * \param rec_names: `string *` Names of the records in the file.
   */
  FileSchema(int num_rec, const std::string *rec_types, const std::string *rec_names=NULL);

  /*! \brief FileSchema instance destroyer.
   */
  ~FileSchema();

  /*! \brief Get the number of records in file.
   *
   * \return num_records: `int` The number of records contained in the file.
   */
  int get_record_number() { return num_records; }

  /*! \brief Get a copy of the record names.
   *
   * \return rec_types: `string *` A new vector of strings with description of records.
   */
  std::string *get_record_names();
  
  /*! \brief Get a copy of the record types.
   *
   * \return rec_names: `string *` A new vector of strings with record names.
   */
  std::string *get_record_types();
};

/*! \class HDFFile
 *
 * \brief HDF5 I/O wrapper class.
 *
 * This class manages I/O operations toward HDF5 format files.
 */
class HDFFile {
 protected:
  //! \brief Identifier list.
  List<hid_t> *id_list;
  //! \brief Name of the file.
  std::string file_name;
  //! \brief Flag for the open file status.
  bool file_open_flag;
  //! File identifier handle.
  hid_t file_id;
  //! Return status of the last operation.
  herr_t status;

 public:
  /*! \brief HDFFile instance constructor.
   *
   * \param name: `string` Name of the file.
   * \param flags: `unsigned int` File access flags (default is `H5F_ACC_EXCL`).
   * \param fcpl_id: `hid_t` File creation property list identifier (default is `H5P_DEFAULT`).
   * \param fapl_id: `hid_t` File access property list identifier (default is `H5P_DEFAULT`).
   */
  HDFFile(
	  const std::string& name, unsigned int flags=H5F_ACC_EXCL,
	  hid_t fcpl_id=H5P_DEFAULT, hid_t fapl_id=H5P_DEFAULT
  );

  /*! \brief HDFFile instance destroyer.
   */
  ~HDFFile();
  
  /*! \brief Close the current file.
   */
  herr_t close();

  /*! \brief Create an empty file from a `FileSchema` instance.
   *
   * \param schema: `FileSchema &` Reference to `FileSchema` instance.
   * \param name: `string` Name of the file.
   * \param flags: `unsigned int` File access flags (default is `H5F_ACC_EXCL`).
   * \param fcpl_id: `hid_t` File creation property list identifier (default is `H5P_DEFAULT`).
   * \param fapl_id: `hid_t` File access property list identifier (default is `H5P_DEFAULT`).
   * \return hdf_file: `HDFFile *` Pointer to a new, open HDF5 file.
   */
  static HDFFile* from_schema(
			      FileSchema &schema, const std::string& name, unsigned int flags=H5F_ACC_EXCL,
			      hid_t fcpl_id=H5P_DEFAULT, hid_t fapl_id=H5P_DEFAULT
  );

  /*! \brief Get current status.
   */
  hid_t get_file_id() { return file_id; }
  
  /*! \brief Get current status.
   */
  herr_t get_status() { return status; }
  
  /*! \brief Check whether the attached file is currently open.
   */
  bool is_open() { return file_open_flag; }
  
  /*! \brief Read data from attached file.
   *
   * \param dataset_name: `string` Name of the dataset to read from.
   * \param data_type: `string` Memory data type identifier.
   * \param buffer: `hid_t` Starting address of the memory sector to store the data.
   * \param mem_space_id: `hid_t` Memory data space identifier (defaults to `H5S_ALL`).
   * \param file_space_id: `hid_t` File space identifier (defaults to `H5S_ALL`).
   * \param dapl_id: `hid_t` Data access property list identifier (defaults to `H5P_DEFAULT`).
   * \param dxpl_id: `hid_t` Data transfer property list identifier (defaults to `H5P_DEFAULT`).
   * \return status: `herr_t` Exit status of the operation.
   */
  herr_t read(
	       const std::string& dataset_name, const std::string& data_type, void *buffer,
	       hid_t mem_space_id=H5S_ALL, hid_t file_space_id=H5S_ALL,
	       hid_t dapl_id=H5P_DEFAULT, hid_t dxpl_id=H5P_DEFAULT
  );

  /*! \brief Write data to attached file.
   *
   * \param dataset_name: `string` Name of the dataset to write to.
   * \param data_type: `string` Memory data type identifier.
   * \param buffer: `hid_t` Starting address of the memory sector to be written.
   * \param mem_space_id: `hid_t` Memory data space identifier (defaults to `H5S_ALL`).
   * \param file_space_id: `hid_t` File space identifier (defaults to `H5S_ALL`).
   * \param dapl_id: `hid_t` Data access property list identifier (defaults to `H5P_DEFAULT`).
   * \param dxpl_id: `hid_t` Data transfer property list identifier (defaults to `H5P_DEFAULT`).
   * \return status: `herr_t` Exit status of the operation.
   */
  herr_t write(
	       const std::string& dataset_name, const std::string& data_type, const void *buffer,
	       hid_t mem_space_id=H5S_ALL, hid_t file_space_id=H5S_ALL,
	       hid_t dapl_id=H5P_DEFAULT, hid_t dxpl_id=H5P_DEFAULT
  );
};

/*! \class VirtualAsciiFile
 *
 * \brief Virtual representation of an ASCII file.
 */
class VirtualAsciiFile {
protected:
  //! \brief The number of lines.
  // int32_t _num_lines;
  //! \brief A vector of strings representing the file lines.
  std::vector<std::string> *_file_lines;

public:
  // const int32_t &num_lines = _num_lines;
  /*! \brief VirtualAsciiFile instance constructor.
   *
   * \param lines: `int32_t` Number of lines, if known in advance (optional, default is 0).
   */
  VirtualAsciiFile(int32_t lines = 0);

  /*! \brief VirtualAsciiFile copy constructor.
   *
   * \param rhs: `const VirtualAsciiFile&` Reference to a VirtualAsciiFile instance.
   */
  VirtualAsciiFile(const VirtualAsciiFile& rhs);

  /*! \brief VirtualAsciiFile instance constructor copying all contents off MPISend() calls from MPI process rr.
   *
   * \param mpidata: `mixMPI *` pointer to MPI data structure.
   * \param rr: `int` rank of the MPI process sending the data.
   */
  VirtualAsciiFile(const mixMPI *mpidata, int rr);

  /*! \brief VirtualAsciiFile instance destroyer.
   */
  ~VirtualAsciiFile();

  /*! \brief Append another VirtualAsciiFile at the end of the current instance.
   *
   * \param rhs: `const VirtualAsciiFile&` Reference to the VirtualAsciiFile to be appended.
   */
  void append(VirtualAsciiFile& rhs);

  /*! \brief Append a line at the end of the file.
   *
   * \param line: `const string&` Reference to a string representing the line.
   */
  void append_line(const std::string& line);

  /*! \brief Append the contents of the VirtualAsciiFile to a physical file on disk.
   *
   * \param file_name: `const string&` Name of the file to append contents to.
   * \return result: `int` A result code (0 if successful).
   */
  int append_to_disk(const std::string& file_name);

  /*! \brief Insert another VirtualAsciiFile at a given position.
   *
   * This function inserts a target VirtualAsciiFile in the current one at the given
   * position. Optionally, a range of lines to be inserted can be specified, otherwise
   * the full content of the target file is inserted. This function DOES NOT increase
   * the size of the inner storage and it can only be used if the inner storage has
   * already been adjusted to contain the insertion target.
   *
   * \param position: `int32_t` The position at which the other file is inserted in this one.
   * \param rhs: `const VirtualAsciiFile&` The refence to the VirtualAsciiFile to be inserted.
   * \param start: `int32_t` The first line to be inserted (optional, default is 0).
   * \param end: `int32_t` The last line to be inserted (optional, default is 0 to read all).
   * \return result: `int` A result code (0 if successful).
   */
  int insert(int32_t position, VirtualAsciiFile& rhs, int32_t start = 0, int32_t end = 0);
  
  /*! \brief Get the number of lines in the current instance.
   *
   * \return size: `int32_t` The number of lines in the VirtualAsciiFile instance.
   */
  int32_t number_of_lines() { return _file_lines->size(); }
    
  /*! \brief Write virtual file contents to a real file on disk.
   *
   * \param file_name: `const string&` Name of the file to append contents to.
   * \return result: `int` A result code (0 if successful).
   */
  int write_to_disk(const std::string& file_name);

  /*! \brief Send VirtualAsciiFile instance to MPI process 0 via MPISend() calls.
   *
   * \param mpidata: `mixMPI *` pointer to MPI data structure.
   */
  void mpisend(const mixMPI *mpidata);
};


/*! \class VirtualBinaryLine
 *
 * \brief Virtual representation of a binary file line.
 */
class VirtualBinaryLine {
// protected:
//   //! \brief The pointer to the piece of data to be written, cast to char *
//   char *_data_pointer;
//   //! \brief the size of the data block.
//   size_t _data_size;

public:
  //! \brief The pointer to the piece of data to be written, cast to char *
  char *_data_pointer;
  //! \brief the size of the data block.
  size_t _data_size;
  //! \brief Read only view of `_data_pointer`.
  const char* data_pointer = _data_pointer;
  //! \brief Read only view of `_data_size`.
  const size_t & data_size = _data_size;

  /*! \brief VirtualBinaryLine instance constructor for `int` data.
   *
   * \param mydata: `int` The piece of data to put in the line.
   */
  VirtualBinaryLine(int mydata);
  
  /*! \brief VirtualBinaryLine instance constructor for `long` data.
   *
   * \param mydata: `long` The piece of data to put in the line.
   */
  VirtualBinaryLine(long mydata);

  /*! \brief VirtualBinaryLine instance constructor for single-precision floating point data.
   *
   * \param mydata: `float` The piece of data to put in the line.
   */
  VirtualBinaryLine(float mydata);

  /*! \brief VirtualBinaryLine instance constructor for double-precision floating point data.
   *
   * \param mydata: `double` The piece of data to put in the line.
   */
  VirtualBinaryLine(double mydata);

  /*! \brief VirtualBinaryLine instance constructor for `dcomplex` data.
   *
   * \param mydata: `dcomplex` The piece of data to put in the line.
   */
  VirtualBinaryLine(dcomplex mydata);

  /*! \brief VirtualBinaryLine copy constructor.
   *
   * \param rhs: `const VirtualBinaryLine&` Reference to a VirtualBinaryLine instance.
   */
  VirtualBinaryLine(const VirtualBinaryLine& rhs);

  /*! \brief VirtualBinaryLine instance constructor copying all contents off MPISend() calls from MPI process rr.
   *
   * \param mpidata: `mixMPI *` pointer to MPI data structure.
   * \param rr: `int` rank of the MPI process sending the data.
   */
  VirtualBinaryLine(const mixMPI *mpidata, int rr);

  /*! \brief VirtualBinaryLine instance destroyer.
   */
  ~VirtualBinaryLine();

  /*! \brief Send VirtualBinaryLine instance to MPI process 0 via MPISend() calls.
   *
   * \param mpidata: `mixMPI *` pointer to MPI data structure.
   */
  void mpisend(const mixMPI *mpidata);
};


/*! \class VirtualBinaryFile
 *
 * \brief Virtual representation of a binary file.
 */
class VirtualBinaryFile {
protected:
  //! \brief The number of lines.
  // int32_t _num_lines;
  // //! \brief A vector of strings representing the file lines.
  // std::vector<VirtualBinaryLine> *_file_lines;

public:
  //! \brief A vector of strings representing the file lines.
  std::vector<VirtualBinaryLine> *_file_lines;
  // const int32_t &num_lines = _num_lines;
  /*! \brief VirtualBinaryFile empty instance constructor.
   *
   */
  VirtualBinaryFile();

  /*! \brief VirtualBinaryFile copy constructor.
   *
   * \param rhs: `const VirtualBinaryFile&` Reference to a VirtualBinaryFile instance.
   */
  VirtualBinaryFile(const VirtualBinaryFile& rhs);

  /*! \brief VirtualBinaryFile instance constructor copying all contents off MPISend() calls from MPI process rr.
   *
   * \param mpidata: `mixMPI *` pointer to MPI data structure.
   * \param rr: `int` rank of the MPI process sending the data.
   */
  VirtualBinaryFile(const mixMPI *mpidata, int rr);

  /*! \brief VirtualBinaryFile instance destroyer.
   */
  ~VirtualBinaryFile();

  /*! \brief Append another VirtualBinaryFile at the end of the current instance.
   *
   * \param rhs: `const VirtualBinaryFile&` Reference to the VirtualBinaryFile to be appended.
   */
  void append(VirtualBinaryFile& rhs);

  /*! \brief Append a line at the end of the file.
   *
   * \param line: `const string&` Reference to a string representing the line.
   */
  void append_line(const VirtualBinaryLine& line);

  /*! \brief Append the contents of the VirtualBinaryFile to a physical file on disk.
   *
   * \param file_name: `const string&` Name of the file to append contents to.
   * \return result: `int` A result code (0 if successful).
   */
  int append_to_disk(const std::string& file_name);

  /*! \brief Get the number of lines in the current instance.
   *
   * \return size: `int32_t` The number of lines in the VirtualBinaryFile instance.
   */
  int32_t number_of_lines() { return _file_lines->size(); }
    
  /*! \brief Write virtual file contents to a real file on disk.
   *
   * \param file_name: `const string&` Name of the file to append contents to.
   * \return result: `int` A result code (0 if successful).
   */
  int write_to_disk(const std::string& file_name);

  /*! \brief Send VirtualBinaryFile instance to MPI process 0 via MPISend() calls.
   *
   * \param mpidata: `mixMPI *` pointer to MPI data structure.
   */
  void mpisend(const mixMPI *mpidata);
};
#endif


