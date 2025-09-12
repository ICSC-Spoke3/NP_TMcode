#include <string>

#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
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

int test_inclusion_hdf5_output();
int test_inclusion_devel();

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
