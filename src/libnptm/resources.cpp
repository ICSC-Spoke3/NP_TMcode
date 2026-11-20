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

/*! \file resources.cpp
 *
 * \brief Implementation of the checks of the memory available to a calculation.
 * \ingroup support
 */

#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#ifndef MPI_VERSION
#include <mpi.h>
#endif
#endif


#ifndef INCLUDE_TYPES_H_
#include "../include/types.h"
#endif

#ifndef INCLUDE_ERRORS_H_
#include "../include/errors.h"
#endif

#ifndef INCLUDE_LOGGING_H_
#include "../include/logging.h"
#endif

#ifndef INCLUDE_LIST_H_
#include "../include/List.h"
#endif

#ifndef INCLUDE_CONFIGURATION_H_
#include "../include/Configuration.h"
#endif

#ifndef INCLUDE_COMMONS_H_
#include "../include/Commons.h"
#endif

#ifndef INCLUDE_RESOURCES_H_
#include "../include/resources.h"
#endif

#ifndef INCLUDE_GPU_RUNTIME_H_
#include "../include/gpu_runtime.h"
#endif

using namespace std;

//! \cond
namespace {
  const double GIB = 1024.0 * 1024.0 * 1024.0;

  // Format a size in bytes as GiB.
  string gib(long bytes) {
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%.3g GiB", bytes / GIB);
    return string(buffer);
  }

  // First value of a one-line system file, or -1 (also for "max").
  long read_value(const string& path) {
    ifstream file(path);
    string token;
    if (!(file >> token)) return -1;
    try {
      size_t used = 0;
      long long value = stoll(token, &used);
      if (used != token.size() || value < 0) return -1;
      return (value > (1LL << 60)) ? -1 : (long)value; // huge values mean no limit
    } catch (const exception &) {
      return -1;
    }
  }

  // Smaller of two sizes, where -1 means unknown.
  long smaller(long a, long b) {
    if (a < 0) return b;
    if (b < 0) return a;
    return (a < b) ? a : b;
  }

  // Free part of the memory limits of the control group of this process and
  // of its parents (cgroup v2 and v1), or -1 if there is no limit.
  long cgroup_available() {
    ifstream file("/proc/self/cgroup");
    string line;
    long result = -1;
    while (getline(file, line)) {
      size_t first = line.find(':');
      size_t second = line.find(':', first + 1);
      if (first == string::npos || second == string::npos) continue;
      string controllers = line.substr(first + 1, second - first - 1);
      string path = line.substr(second + 1);
      string root, limit_file, usage_file;
      if (line.compare(0, 3, "0::") == 0) {
	root = "/sys/fs/cgroup";
	limit_file = "/memory.max";
	usage_file = "/memory.current";
      } else if (("," + controllers + ",").find(",memory,") != string::npos) {
	root = "/sys/fs/cgroup/memory";
	limit_file = "/memory.limit_in_bytes";
	usage_file = "/memory.usage_in_bytes";
      } else {
	continue;
      }
      // Walk from the group of the process up to the root of the hierarchy.
      while (true) {
	string folder = root + ((path == "/") ? "" : path);
	long limit = read_value(folder + limit_file);
	long usage = read_value(folder + usage_file);
	if (limit >= 0 && usage >= 0) {
	  result = smaller(result, (limit > usage) ? limit - usage : 0);
	}
	if (path.empty() || path == "/") break;
	size_t slash = path.find_last_of('/');
	path = (slash == 0 || slash == string::npos) ? "/" : path.substr(0, slash);
      }
    }
    return result;
  }

  // MemAvailable of /proc/meminfo, or -1.
  long meminfo_available() {
    ifstream file("/proc/meminfo");
    string line;
    while (getline(file, line)) {
      if (line.compare(0, 13, "MemAvailable:") == 0) {
	istringstream fields(line.substr(13));
	long kib;
	if (fields >> kib) return kib * 1024;
      }
    }
    return -1;
  }

  // Number of OpenMP threads of the outer parallel regions.
  int count_threads() {
    int threads = 1;
#ifdef _OPENMP
#pragma omp parallel
    {
      if (omp_get_thread_num() == 0) threads = omp_get_num_threads();
    }
#endif
    return threads;
  }

  // Advice for a calculation that does not fit.
  const char *PARALLEL_HINT = "reduce the model size or the number of OpenMP threads or MPI "
    "processes per node, or set MEMORY_POLICY=FORCE_LARGE_RUN in the input to run anyway";
  const char *SERIAL_HINT = "reduce the model size, or set MEMORY_POLICY=FORCE_LARGE_RUN "
    "in the input to run anyway";

  // Stop or warn about a shortage, according to the policy.
  void report_shortage(const string& problem, short policy, Logger *logger, const char *hint) {
    if (policy == MEMORY_POLICY_FORCE_LARGE_RUN) {
      if (logger != NULL) {
	logger->log(problem + "; running anyway (MEMORY_POLICY=FORCE_LARGE_RUN)\n", LOG_WARN);
      }
    } else {
      throw MemoryError(problem + ": " + hint);
    }
  }
}
//! \endcond

short parse_memory_policy(const std::string& value) {
  string word = value;
  while (!word.empty() && (word.back() == ' ' || word.back() == '\t' || word.back() == '\r')) word.pop_back();
  if (word == "STOP") return MEMORY_POLICY_STOP;
  if (word == "FORCE_LARGE_RUN") return MEMORY_POLICY_FORCE_LARGE_RUN;
  throw ParseError("unknown memory policy \"" + word + "\" (use STOP or FORCE_LARGE_RUN)");
}

short parse_device_policy(const std::string& value) {
  string word = value;
  while (!word.empty() && (word.back() == ' ' || word.back() == '\t' || word.back() == '\r')) word.pop_back();
  if (word == "STOP") return DEVICE_POLICY_STOP;
  if (word == "USE_CPU_FALLBACK") return DEVICE_POLICY_USE_CPU_FALLBACK;
  throw ParseError("unknown device policy \"" + word + "\" (use STOP or USE_CPU_FALLBACK)");
}

long get_available_host_memory() {
  long result = smaller(meminfo_available(), cgroup_available());
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE)
  if (result < 0) {
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);
    if (pages > 0 && page_size > 0) result = pages * page_size;
  }
#endif
  return result;
}

long get_available_device_memory(int device) {
  return gpu_free_memory(device);
}

long MemoryDemand::host_bytes(bool main_process, int threads) const {
  long result = threads * per_thread;
  if (main_process) {
    // The results of all the scales, and one scale for each other thread.
    result += full_output + (threads - 1) * scale_output;
  } else {
    // Thread 0 collects the scales of all the threads, which hold one each.
    result += (2 * threads - 1) * scale_output;
  }
  return result;
}

ResourceSurvey::ResourceSurvey(const mixMPI *mpidata, int device_count, int status) {
#if defined(USE_TARGET_OFFLOAD) && defined(_OPENMP)
  // Without MAGMA or cuBLAS, the GPUs are the OpenMP offload devices.
  if (device_count == 0) device_count = omp_get_num_devices();
#endif
  const int fields = 6;
  long local[fields];
  local[0] = mpidata->rank;
  local[1] = count_threads();
  local[2] = (device_count > 0) ? mpidata->rank % device_count : -1;
  local[3] = get_available_host_memory();
  local[4] = (local[2] >= 0) ? get_available_device_memory((int)local[2]) : -1;
  local[5] = status;
  int nprocs = 1;
  vector<long> all(local, local + fields);
#ifdef MPI_VERSION
  if (mpidata->mpirunning) {
    // The processes that share memory belong to one node: name it by the
    // lowest rank among them.
    MPI_Comm node_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpidata->rank, MPI_INFO_NULL, &node_comm);
    int leader = mpidata->rank;
    MPI_Bcast(&leader, 1, MPI_INT, 0, node_comm);
    MPI_Comm_free(&node_comm);
    local[0] = leader;
    nprocs = mpidata->nprocs;
    all.assign((size_t)fields * nprocs, 0);
    MPI_Gather(local, fields, MPI_LONG, all.data(), fields, MPI_LONG, 0, MPI_COMM_WORLD);
  }
#endif
  for (int p = 0; p < nprocs; p++) {
    node.push_back((int)all[fields * p]);
    threads.push_back((int)all[fields * p + 1]);
    device.push_back((int)all[fields * p + 2]);
    host_free.push_back(all[fields * p + 3]);
    device_free.push_back(all[fields * p + 4]);
    device_status.push_back((int)all[fields * p + 5]);
  }
}

void check_memory(
  const ResourceSurvey& survey, const MemoryDemand& demand, double host_ram_gb,
  double gpu_ram_gb, short policy, Logger *logger
) {
  // Needs and availability by node and by (node, GPU).
  map<int, long> node_need, node_free;
  map<pair<int, int>, long> gpu_need, gpu_free;
  map<pair<int, int>, int> gpu_processes;
  for (size_t p = 0; p < survey.node.size(); p++) {
    int node = survey.node[p];
    node_need[node] += demand.host_bytes(p == 0, survey.threads[p]);
    node_free[node] = smaller(node_free.count(node) ? node_free[node] : -1, survey.host_free[p]);
    if (survey.device[p] >= 0 && demand.device_per_thread > 0) {
      pair<int, int> gpu(node, survey.device[p]);
      gpu_need[gpu] += survey.threads[p] * demand.device_per_thread + NPTM_GPU_PROCESS_OVERHEAD;
      gpu_processes[gpu]++;
      gpu_free[gpu] = smaller(gpu_free.count(gpu) ? gpu_free[gpu] : -1, survey.device_free[p]);
    }
  }
  const bool multi_node = (node_need.size() > 1);
  for (const auto &item : node_need) {
    string where = multi_node ? " on the node of MPI process " + to_string(item.first) : "";
    long available;
    string source;
    if (host_ram_gb > 0.0) {
      available = (long)(host_ram_gb * GIB);
      source = "allowed by HOST_RAM_GB";
    } else {
      available = node_free[item.first];
      source = "available";
    }
    if (available < 0) {
      if (logger != NULL) {
	logger->log(
	  "the available host memory" + where + " is unknown: the calculation needs "
	  + gib(item.second) + " and its size is not checked (set HOST_RAM_GB to check it).\n",
	  LOG_WARN
	);
      }
      continue;
    }
    if (logger != NULL) {
      logger->log(
	"the calculation needs " + gib(item.second) + " of host memory" + where + ", "
	+ gib(available) + " " + source + ".\n", LOG_INFO
      );
    }
    if (item.second > available) {
      report_shortage(
	"the calculation needs " + gib(item.second) + " of host memory" + where
	+ ", more than the " + gib(available) + " " + source, policy, logger, PARALLEL_HINT
      );
    }
  }
  for (const auto &item : gpu_need) {
    string where = "GPU " + to_string(item.first.second)
      + (multi_node ? " of the node of MPI process " + to_string(item.first.first) : "");
    long available;
    string source;
    if (gpu_ram_gb > 0.0) {
      available = (long)(gpu_ram_gb * GIB);
      source = "allowed by GPU_RAM_GB";
    } else if (gpu_free[item.first] >= 0) {
      // A margin for the allocations that the estimates do not cover.
      available = (long)(NPTM_GPU_USABLE_FRACTION * gpu_free[item.first]);
      source = "usable (" + to_string((int)(100 * NPTM_GPU_USABLE_FRACTION + 0.5)) + "% of "
	+ gib(gpu_free[item.first]) + " free)";
    } else {
      available = -1;
      source = "free";
    }
    if (available < 0) {
      if (logger != NULL) {
	logger->log(
	  "the free memory of " + where + " is unknown: the calculation needs "
	  + gib(item.second) + " and its size is not checked (set GPU_RAM_GB to check it).\n",
	  LOG_WARN
	);
      }
      continue;
    }
    if (logger != NULL) {
      logger->log(
	"the calculation needs " + gib(item.second) + " of memory on " + where + " (with "
	+ gib(gpu_processes[item.first] * NPTM_GPU_PROCESS_OVERHEAD) + " for the GPU runtime).\n", LOG_INFO
      );
      logger->log(gib(available) + " " + source + ".\n", LOG_INFO);
    }
    if (item.second > available) {
      report_shortage(
	"the calculation needs " + gib(item.second) + " of memory on " + where
	+ ", more than the " + gib(available) + " " + source, policy, logger, PARALLEL_HINT
      );
    }
  }
}

void check_devices(const ResourceSurvey& survey, short policy, Logger *logger) {
  string problems;
  int count = 0;
  for (size_t p = 0; p < survey.device_status.size(); p++) {
    int status = survey.device_status[p];
    if (status == DEVICE_STATUS_OK) continue;
    string what = (status == DEVICE_STATUS_NONE_FOUND) ? "found no GPU" : "could not initialize or use its GPU";
    if (count > 0) problems += "; ";
    problems += ((survey.device_status.size() > 1) ? "MPI process " + to_string(p) + " " : "the program ") + what;
    count++;
  }
  if (count == 0) return;
  if (policy == DEVICE_POLICY_USE_CPU_FALLBACK) {
    if (logger != NULL) {
      logger->log(problems + ": inverting the matrices on the CPU (DEVICE_POLICY=USE_CPU_FALLBACK).\n", LOG_WARN);
    }
  } else {
    throw DeviceError(
      problems + ": check the GPUs and the CUDA installation, or set "
      "DEVICE_POLICY=USE_CPU_FALLBACK in the input to use the CPU"
    );
  }
}

void check_host_memory(long bytes, const std::string& what, double host_ram_gb, short policy, Logger *logger) {
  long available;
  string source;
  if (host_ram_gb > 0.0) {
    available = (long)(host_ram_gb * GIB);
    source = "allowed by HOST_RAM_GB";
  } else {
    available = get_available_host_memory();
    source = "available";
  }
  if (available < 0) {
    if (logger != NULL) {
      logger->log(
	"the available host memory is unknown: " + what + " need " + gib(bytes)
	+ " and their size is not checked (set HOST_RAM_GB to check it).\n", LOG_WARN
      );
    }
    return;
  }
  if (logger != NULL) {
    logger->log(what + " need " + gib(bytes) + " of host memory, " + gib(available) + " " + source + ".\n", LOG_INFO);
  }
  if (bytes > available) {
    report_shortage(
      what + " need " + gib(bytes) + " of host memory, more than the " + gib(available) + " " + source,
      policy, logger, SERIAL_HINT
    );
  }
}
