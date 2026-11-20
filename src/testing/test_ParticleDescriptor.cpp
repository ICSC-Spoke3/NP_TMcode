/*! \file test_ParticleDescriptor.cpp
 *
 * \brief Test the construction of `ParticleDescriptor` instances.
 *
 * Builds a sphere, a cluster and an inclusion particle descriptor from
 * the legacy test data in `test_data` and checks that no exception is
 * raised. Must be run from a folder two levels below the repository
 * root (e.g. `build/testing`), since data paths are relative.
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

using namespace std;

/*! \brief Test a `ParticleDescriptorCluster` built from `test_data/cluster/case_3`.
 *
 * Read the legacy configuration files of the test case and build
 * the corresponding particle descriptor, then release it.
 *
 * \return result: `int` 0 if the test passed, 100 if it failed.
 */
int test_cluster_case_3();

/*! \brief Test a `ParticleDescriptorCluster` built from `test_data/cluster`.
 *
 * Read the legacy configuration files of the test case and build
 * the corresponding particle descriptor, then release it.
 *
 * \return result: `int` 0 if the test passed, 10 if it failed.
 */
int test_cluster_devel();

/*! \brief Test a `ParticleDescriptorInclusion` built from `test_data/inclusion`.
 *
 * Read the legacy configuration files of the test case and build
 * the corresponding particle descriptor, then release it.
 *
 * \return result: `int` 0 if the test passed, 1000 if it failed.
 */
int test_inclusion();

/*! \brief Test a `ParticleDescriptorSphere` built from `test_data/sphere`.
 *
 * Read the legacy configuration files of the test case and build
 * the corresponding particle descriptor, then release it.
 *
 * \return result: `int` 0 if the test passed, 1 if it failed.
 */
int test_sphere();

/*! \brief Main program execution body.
 *
 * Runs each test in sequence and adds up their failure codes, so that
 * every decimal digit of the result identifies one failing test.
 *
 * \return result: `int` Sum of the failure codes (1, 10, 100 and 1000); 0 if all tests passed.
 */
int main() {
  int result = 0;
  result += test_sphere(); // 1 if failed
  result += test_cluster_devel(); // 10 if failed
  result += test_cluster_case_3(); // 100 if failed
  result += test_inclusion(); // 1000 if failed
  return result;
}

int test_cluster_case_3() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/cluster/case_3/DCLU";
    const string scat_data_file = "../../test_data/cluster/case_3/DEDFB_33";
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    ParticleDescriptor *pd = new ParticleDescriptorCluster(gconf, sconf);
    delete gconf;
    delete sconf;
    delete pd;
  } catch (const exception& ex) {
    result = 100;
  }
  return result;
}

int test_cluster_devel() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/cluster/DCLU";
    const string scat_data_file = "../../test_data/cluster/DEDFB";
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    ParticleDescriptorCluster *pd = new ParticleDescriptorCluster(gconf, sconf);
    delete gconf;
    delete sconf;
    delete pd;
  } catch (const exception& ex) {
    result = 10;
  }
  return result;
}

int test_inclusion() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/inclusion/DINCLU";
    const string scat_data_file = "../../test_data/inclusion/DEDFB";
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    ParticleDescriptorInclusion *pd = new ParticleDescriptorInclusion(gconf, sconf);
    delete gconf;
    delete sconf;
    delete pd;
  } catch (const exception& ex) {
    result = 1000;
  }
  return result;
}

int test_sphere() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/sphere/DSPH";
    const string scat_data_file = "../../test_data/sphere/DEDFB";
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    ParticleDescriptorSphere *pd = new ParticleDescriptorSphere(gconf, sconf);
    delete gconf;
    delete sconf;
    delete pd;
  } catch (const exception& ex) {
    result = 1;
  }
  return result;
}
