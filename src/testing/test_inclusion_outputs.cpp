/*! \file test_inclusion_outputs.cpp
 *
 * \brief Test reading and writing of `InclusionOutputInfo` instances.
 *
 * Must be run from a folder two levels below the repository root
 * (e.g. `build/testing`), since data paths are relative.
 * \ingroup tests
 */

#include <string>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

#ifndef INCLUDE_OUTPUTS_H_
#include "../include/outputs.h"
#endif

using namespace std;

/*! \brief Convert `test_data/inclusion/c_OINCLU.hd5` to a legacy `c_OINCLU` file.
 *
 * Load the output of the test case from its HDF5 file and write it
 * again in the legacy text format.
 *
 * \return result: `int` 0 if the test passed, 1 if it failed.
 */
int test_inclusion_hdf5_output();

/*! \brief Allocate an `InclusionOutputInfo` for the `test_data/inclusion` case.
 *
 * Read the legacy configuration files of the test case and allocate
 * an empty output structure for it.
 *
 * \return result: `int` 0 if the test passed, 10 if it failed.
 */
int test_inclusion_devel();

/*! \brief Main program execution body.
 *
 * Runs each test in sequence and adds up their failure codes, so that
 * every decimal digit of the result identifies one failing test.
 *
 * \return result: `int` Sum of the failure codes (1 and 10); 0 if all tests passed.
 */
int main() {
  int result = 0;
  result += test_inclusion_hdf5_output(); // 1 if failed
  result += test_inclusion_devel(); // 10 if failed
  return result;
}

int test_inclusion_hdf5_output() {
  int result = 0;
  try {
    const string hdf5_file = "../../test_data/inclusion/c_OINCLU.hd5";
    InclusionOutputInfo *oi = new InclusionOutputInfo(hdf5_file);
    oi->write("c_OINCLU", "LEGACY");
    delete oi;
  } catch (const exception& ex) {
    result = 1;
  }
  return result;
}

int test_inclusion_devel() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/inclusion/DINCLU";
    const string scat_data_file = "../../test_data/inclusion/DEDFB";
    mixMPI *mpidata = new mixMPI();
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    InclusionOutputInfo *oi = new InclusionOutputInfo(sconf, gconf, mpidata);
    delete gconf;
    delete sconf;
    delete oi;
    delete mpidata;
  } catch (const exception& ex) {
    result = 10;
  }
  return result;
}
