/*! \file test_sphere_outputs.cpp
 *
 * \brief Test reading and writing of `SphereOutputInfo` instances.
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

/*! \brief Convert `test_data/sphere/c_OSPH.hd5` to a legacy `c_OSPH` file.
 *
 * Load the output of the test case from its HDF5 file and write it
 * again in the legacy text format.
 *
 * \return result: `int` 0 if the test passed, 1 if it failed.
 */
int test_sphere_hdf5_output();

/*! \brief Allocate a `SphereOutputInfo` for the `test_data/sphere` case.
 *
 * Read the legacy configuration files of the test case and allocate
 * an empty output structure for it.
 *
 * \return result: `int` 0 if the test passed, 10 if it failed.
 */
int test_sphere_devel();

/*! \brief Main program execution body.
 *
 * Runs each test in sequence and adds up their failure codes, so that
 * every decimal digit of the result identifies one failing test.
 *
 * \return result: `int` Sum of the failure codes (1 and 10); 0 if all tests passed.
 */
int main() {
  int result = 0;
  result += test_sphere_hdf5_output(); // 1 if failed
  result += test_sphere_devel(); // 10 if failed
  return result;
}

int test_sphere_hdf5_output() {
  int result = 0;
  try {
    const string hdf5_file = "../../test_data/sphere/c_OSPH.hd5";
    SphereOutputInfo *oi = new SphereOutputInfo(hdf5_file);
    oi->write("c_OSPH", "LEGACY");
    delete oi;
  } catch (const exception& ex) {
    result = 1;
  }
  return result;
}

int test_sphere_devel() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/sphere/DSPH";
    const string scat_data_file = "../../test_data/sphere/DEDFB";
    mixMPI *mpidata = new mixMPI();
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    SphereOutputInfo *oi = new SphereOutputInfo(sconf, gconf, mpidata);
    delete gconf;
    delete sconf;
    delete oi;
    delete mpidata;
  } catch (const exception& ex) {
    result = 10;
  }
  return result;
}
