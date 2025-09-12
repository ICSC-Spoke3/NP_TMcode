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

int test_cluster_hdf5_output();
int test_cluster_case_3();

int main() {
  int result = 0;
  result += test_cluster_hdf5_output(); // 1 if failed
  // result += test_cluster_devel(); // 10 if failed
  result += test_cluster_case_3(); // 100 if failed
  // result += test_inclusion(); // 1000 if failed
  return result;
}

int test_cluster_hdf5_output() {
  int result = 0;
  try {
    const string hdf5_file = "../../test_data/cluster/c_OCLU_24.hd5";
    ClusterOutputInfo *oi = new ClusterOutputInfo(hdf5_file);
    oi->write("c_OCLU_24", "LEGACY");
    delete oi;
  } catch (const exception& ex) {
    result = 1;
  }
  return result;
}

int test_cluster_case_3() {
  int result = 0;
  try {
    const string geom_data_file = "../../test_data/cluster/case_3/DCLU";
    const string scat_data_file = "../../test_data/cluster/case_3/DEDFB_33";
    mixMPI *mpidata = new mixMPI();
    GeometryConfiguration *gconf = GeometryConfiguration::from_legacy(geom_data_file);
    ScattererConfiguration *sconf = ScattererConfiguration::from_dedfb(scat_data_file);
    ClusterOutputInfo *oi = new ClusterOutputInfo(sconf, gconf, mpidata);
    delete gconf;
    delete sconf;
    delete oi;
    delete mpidata;
  } catch (const exception& ex) {
    result = 100;
  }
  return result;
}
