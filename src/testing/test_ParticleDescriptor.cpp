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

using namespace std;

int test_cluster_case_3();
int test_cluster_devel();
int test_inclusion();
int test_sphere();

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
