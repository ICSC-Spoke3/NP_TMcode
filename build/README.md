# Folder instructions
This directory contains the configuration and building system.

## Instructions
The set-up and building process of *NP_TMcode* is based on tha automatic detection of the system capabilities, through an *autoconf* generated `configure` script. In case the user wants to override the default configuration, several options are available to modify the configuration behaviour. Use:

   > ./configure --help

to obtain a list of the available configuration options. We refer to the code inline documentation and to the PDF release notes of the distribution for advanced instructions on code set up and execution. In the following sections, we give instructions for a simple installation and some elementary tests.

Use of this code and of any derived implementation should credit the original authors. Scientific publications should do so by citing the following references:

- Saija et al. 2001, ApJ, 559, 993, DOI:10.1086/322350
- Borghese, Denti, Saija 2007, Scattering from Model Nonspherical Particles (ISBN 978-3-540-37413-8), DOI:10.1007/978-3-540-37414-5
- La Mura et al. 2026, Astronomy and Computing, 55, DOI:10.1016/j.ascom.2026.101077

## Installation
The proper way to install the code and execute it is based on pre-configured builds handled by *GNU make*. To install the code on your local machine, you should enter the `build` folder and issue:

   > ./configure [CONFIGURE_OPTIONS]
   >
   > make

where `CONFIGURE_OPTIONS` are optional arguments used to enable and disable specific code features, such as parallel implementation and GPU offload of calculations. These features are enabled by default, if the local system has the hardware and software ability to support them. See the help pages of the `configure` script for more details on the available optional features.

*NOTE on GPU offload:* the OpenMP offload of `np_trapping` can be compiled with `g++` or with `clang++`. With `g++`, each OpenMP thread of a GPU team runs on a whole warp, so the offloaded loops use only a small part of the GPU; `clang++` uses all the threads of the warps and is the recommended compiler for offload builds. To build with `clang++` through an MPI wrapper, select it with the wrapper's environment variable before running `configure` (`OMPI_CXX` for Open MPI, `MPICH_CXX` for MPICH); `configure` records it in `Makefile.in`, so `make` uses the same compiler:

   > export OMPI_CXX=clang++
   > ./configure
   > make

Without MPI, use `CXX=clang++ ./configure`. With `clang++`, `configure` detects the GPU architecture with `nvidia-smi` and looks for a CUDA installation (needed for `ptxas`); use `--with-offload-arch=ARCH` (e.g. `sm_89`) and `--with-offload-cuda=PATH` to set them, for instance when building on a machine without a GPU. With `g++`, the default architecture is `sm_70`.

*NOTE on cuBLAS builds:* without MAGMA, the matrices are inverted with cuSOLVER, which comes with the CUDA toolkit (version 11.1 or later). `configure` enables cuBLAS only if it also finds cuSOLVER.

*NOTE on AMD GPUs (experimental):* MAGMA built on HIP/ROCm and OpenMP offload to AMD GPUs (`--with-offload-target=amdgcn`) are supported by `configure` and by the code, but they have not yet been tested on AMD hardware: only the NVIDIA (CUDA) builds are tested. `configure` reads the MAGMA back end from `magma_config.h` and uses the ROCm installation in `ROCM_PATH` (default `/opt/rocm`); the AMD architecture (e.g. `gfx90a`) is detected with `amdgpu-arch` or `rocminfo`, or set with `--with-offload-arch=ARCH`. cuBLAS is not available with a MAGMA built on HIP. Please report any problem found with these builds.

## Code version

The version of the package is the `version` key of `CITATION.cff`, in the root folder of the repository. `configure` reads it ("checking the code version...") and passes it to the programs, which print it with `--version`; after changing it, run `configure` again. The Python scripts and the user's manual read the same key.

## Code work-flow
The original code produces output in the current working directory (the path where the code is executed from). The build directory is intended to collect local builds and test run output in a safe place, without cluttering the code development folders, thus keeping the source folders as clean as possible.

This section describes the use of the single programs that make up the *NP_TMcode* suite, once the binaries have been properly built.

*NOTE on MPI builds:* `configure` enables MPI when the C++ compiler is an MPI wrapper, such as `mpicxx` (the build then has `-DUSE_MPI` in `Makefile.in`). In an MPI build, `np_sphere`, `np_cluster` and `np_inclusion` must be started through the MPI launcher, even to run a single process:

   > mpirun -n 1 ./np_cluster

Started directly, they may hang at start-up, depending on the MPI library. Use `-n N` to run on `N` processes; the number of OpenMP threads of each process is set with the `OMP_NUM_THREADS` environment variable. `np_trapping` does not use MPI and can always be started directly. The commands below are written for a build without MPI. Please, refer to the code user manual in the `doc` folder for more details on parallel runs.

### cluster

1. cd to the `build/cluster` folder.
2. copy the test data in the current working folder

   > cp ../../test_data/cluster/DEDFB .
   >
   > cp ../../test_data/cluster/DCLU .
   
3. Run `edfb_clu`:
   
   > ./edfb_clu
   
4. Run `clu`:
   
   > ./clu
   
*NOTE:* both `edfb_clu` and `clu` expect to find their input in the current working directory. The input used in this example is located in a folder named `../../test_data/cluster/` (i.e. two levels above the current execution path).

5. Run `np_cluster` (`mpirun -n 1 ./np_cluster` in an MPI build):

   > ./np_cluster

*NOTE:* The *C++* version does not need to run a configuration program because all configuration operations are handled by the code at run-time.

6. Check the consistency between the output files (the default output file for the *FORTRAN* code is named `OCLU`, while the corresponding *C++* output has the default name of `c_OCLU`). This can be done, e. g., using the `pycompare.py` script (which requires _Python 3_):

   > ../../src/scripts/pycompare.py --ffile OCLU --cfile c_OCLU

The default behaviour of `np_cluster` is to take the same input files as `edfb_clu` and `clu` and to write the output in the current folder. If needed, different input and output paths can be given as command-line arguments:

   > ./np_cluster PATH_TO_DEDFB PATH_TO_DCLU OUTPUT_PATH

### sphere

1. cd to the `build/sphere` folder.
2. copy the test data in the current working folder

   > cp ../../test_data/sphere/DEDFB .
   >
   > cp ../../test_data/sphere/DSPH .
   
3. Run `edfb_sph`:

   > ./edfb_sph
   
4. Run `sph`:

   > ./sph
   
*NOTE:* both `edfb_sph` and `sph` expect to find their input in the current working directory. The input used in this example is located in a folder named `../../test_data/sphere/` (i.e. two levels above the current execution path).

5. Run `np_sphere` (`mpirun -n 1 ./np_sphere` in an MPI build):

   > ./np_sphere

*NOTE:* The *C++* version does not need to run a configuration program because all configuration operations are handled by the code at run-time.

5. Check the consistency between the output files (the default output file for the *FORTRAN* code is named `OSPH`, while the corresponding *C++* output has the default name of `c_OSPH`). This can be done, e. g., using the `pycompare.py` script (which requires _Python 3_):

   > ../../src/scripts/pycompare.py --ffile OSPH --cfile c_OSPH

The default behaviour of `np_sphere` is to take the same input files as `edfb_sph` and `sph` and to write the output in the current folder. If needed, different input and output paths can be given as command-line arguments:

   > ./np_sphere PATH_TO_DEDFB PATH_TO_DSPH OUTPUT_PATH

### trapping

The execution of trapping programs requires at least one of the previous programs to have produced a complete output set. A light-weight trapping calculation has been configured with input and legacy output files stored in the `../../test_data/trapping/` folder. Since the FORTRAN code assumes the input and output to be defined within the program, it is not yet possible to run the *FORTRAN* version on this case, unless the source code is modified accordingly. Conversely, the *C++* version can be executed without the need to modify and re-compile the code. The work-flow to test trapping is described below.

1. cd to the `build/sphere` folder.
2. run `np_sphere` with arguments to take input from the trapping test and write output in the trapping build folder (in an MPI build, start it as `mpirun -n 1 ./np_sphere` followed by the same arguments):

   > ./np_sphere ../../test_data/trapping/DEDFB ../../test_data/trapping/DSPH ../trapping

3. cd to the trapping folder.
4. run `np_trapping`

   > ./np_trapping ../../test_data/trapping/DFRFME ../../test_data/trapping/DLFFFT .

5. Check the consistency between `np_trapping` output files (`c_force_cs.txt` and `c_torque_cs.txt`) and the legacy *FORTRAN* output for this case (the files named, respectively, `fort.66` and `fort.67` in the `test_data/trapping/` folder). Consider that some of the output values will be affected by numeric noise and take substantially different values. However, this is expected for results whose order of magnitude is clearly below the precision level of the calculation, as they represent results approaching zero that were just approximated with different precision. This can be done using the `pycompare.py` script (which requires _Python 3_):

   > ../../src/scripts/pycompare.py --ffile ../../test_data/trapping/fort.66 --cfile c_force_cs.txt
   >
   > ../../src/scripts/pycompare.py --ffile ../../test_data/trapping/fort.67 --cfile c_torque_cs.txt --data-order=-10

(where the `--data-order=-10` option tells the script to consider as data only nummbers down to the order of 1e-10 and ignore anything below this limit).

### testing

The `testing` folder contains programs that are used to test the consistency of data written in different binary formats. At present, the only binary files that can be read and written both in legacy and _HDF5_ format are the configuration files (named `TEDF`) and the transition matrix files (named `TTMS`). The _HDF5_ versions of these files (marked with an `hd5` extension) can be inspected by standard _HDF5_ tools, such as `h5dump`. To check the consistency of these two types of files, the `testing` suite provides two programs:

- `test_TEDF`

This program checks for the consistency of the configuration data loaded from a formatted configuration file, a legacy binary file and an _HDF5_ binary file. It must be executed as:

    > ./test_TEDF PATH_TO_EDFB PATH_TO_c_TEDF PATH_TO_c_TEDF.hd5

where the command line arguments must be valid paths to, respectively, the formatted configuration file, a binary configuration file produced by one of the _NP_ codes presented above (generally named `c_TEDF`) and a binary configuration file saved in _HDF5_ format (`c_TEDF.hd5`). The program checks for the data in each of the above formats and, after writing a short log to the terminal, returns a result value of 0, in case of success, or a positive error code in case of inconsistent data.

- `test_TTMS`

This program checks for the consistency of transition matrix data files. It is executed similarly to the previous one, but with just two arguments:

    > ./test_TTMS PATH_TO_c_TTMS PATH_TO_c_TTMS.hd5

where the arguments must be valid paths to binary transition matrix files saved in the legacy or the _HDF5_ format. Similarly to the previous case, the program writes a short log and returns a result value of 0, for success, or a positive number in case of inconsistency detection.

The `test_memory_estimates` program checks that the memory estimates of the data classes (their `get_size()` methods), which the programs use to report and check their memory requirements, agree with the memory that the classes actually allocate. From the `build/testing` folder, run:

    > ./test_memory_estimates [PATH_TO_test_data]

The path defaults to `../../test_data`. The program builds the sphere, cluster, inclusion and trapping data structures from the test input files, prints the measured and estimated sizes and returns 0 if every estimate is within 1% (or 16 KiB) of the measured size, or 1 otherwise. In MPI builds, start it with `mpirun -n 1`.

# License

   Copyright (C) 2025   INAF - Osservatorio Astronomico di Cagliari

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
