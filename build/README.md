# Folder instructions
This directory contains the configuration and building system.

## Instructions
The set-up and building process of *NP_TMcode* is based on tha automatic detection of the system capabilities, through an *autoconf* generated `configure` script. In case the user wants to override the default configuration, several options are available to modify the configuration behaviour. Use:

   > ./configure --help

to obtain a list of the available configuration options. We refer to the code inline documentation and to the PDF release notes of the distribution for advanced instructions on code set up and execution. In the following sections, we give instructions for a simple installation and some elementary tests.

Use of this code and of any derived implementation should credit the original authors. Scientific publications should do so by citing the following references:

- Saija et al. 2001, ApJ, 559, 993, DOI:10.1086/322350
- Borghese, Denti, Saija 2007, Scattering from Model Nonspherical Particles (ISBN 978-3-540-37413-8), DOI:10.1007/978-3-540-37414-5

## Installation
The proper way to install the code and execute it is based on pre-configured builds handled by *GNU make*. To install the code on your local machine, you should enter the `build` folder and issue:

   > ./configure [CONFIGURE_OPTIONS]

where `CONFIGURE_OPTIONS` are optional arguments used to enable and disable specific code features, such as parallel implementation and GPU offload of calculations. These features are enabled by default, if the local system has the hardware and software ability to support them. See the help pages of the `configure` script for more details on the available optional features.

## Code work-flow
The original code produces output in the current working directory (the path where the code is executed from). The build directory is intended to collect local builds and test run output in a safe place, without cluttering the code development folders, thus keeping the source folders as clean as possible.

This section describes the use of the single programs that make up the *NP_TMcode* suite, once the binaries have been properly built.

*NOTE:* This set of instructions applies to the serial code implementation only. Please, refer to the Release Notes of `NP_TMcode-M8.03` for more details on how to build and use a parallel implementation of the code.

### cluster

1. cd to the `build/cluster` folder.
2. Run `edfb_clu`:
   
   > ./edfb_clu
   
3. Run `clu`:
   
   > ./clu
   
*NOTE:* both `edfb_clu` and `clu` expect an input which is assumed to be in a folder named `../../test_data/cluster/` (i.e. two levels above the current execution path)

4. Run `np_cluster`:

   > ./np_cluster

*NOTE:* The *C++* version does not need to run a configuration program because all configuration operations are handled by the code at run-time.

5. Check the consistency between the output files (the default output file for the *FORTRAN* code is named `OCLU`, while the corresponding *C++* output has the default name of `c_OCLU`).

The default behaviour of `np_cluster` is to take the same input files as `edfb_clu` and `clu` and to write the output in the current folder. If needed, different input and output paths can be given as command-line arguments:

   > ./np_cluster PATH_TO_DEDFB PATH_TO_DCLU OUTPUT_PATH

### sphere

1. cd to the `build/sphere` folder.
2. Run `edfb_sph`:

   > ./edfb_sph
   
3. Run `sph`:

   > ./sph
   
*NOTE:* both `edfb_sph` and `sph` expect an input which is assumed to be in a folder named `../../test_data/sphere/` (i.e. two levels above the current execution path)

4. Run `np_sphere`:

   > ./np_sphere

*NOTE:* The *C++* version does not need to run a configuration program because all configuration operations are handled by the code at run-time.

5. Check the consistency between the output files (the default output file for the *FORTRAN* code is named `OSPH`, while the corresponding *C++* output has the default name of `c_OSPH`).

The default behaviour of `np_sphere` is to take the same input files as `edfb_sph` and `sph` and to write the output in the current folder. If needed, different input and output paths can be given as command-line arguments:

   > ./np_sphere PATH_TO_DEDFB PATH_TO_DSPH OUTPUT_PATH

### trapping

The execution of trapping programs requires at least one of the previous programs to have produced a complete output set. A light-weight trapping calculation has been configured with input and legacy output files stored in the `../../test_data/trapping/` folder. Since the FORTRAN code assumes the input and output to be defined within the program, it is not yet possible to run the *FORTRAN* version on this case, unless the source code is modified accordingly. Conversely, the *C++* version can be executed without the need to modify and re-compile the code. The work-flow to test trapping is described below.

1. cd to the `build/sphere` folder.
2. run `np_sphere` with arguments to take input from the trapping test and write output in the trapping build folder:

   > ./np_sphere ../../test_data/trapping/DEDFB ../../test_data/trapping/DSPH ../trapping

3. cd to the trapping folder.
4. run `np_trapping`

   > ./np_trapping ../../test_data/trapping/DFRFME ../../test_data/trapping/DLFFT .

5. Check the consistency between `np_trapping` output files (`c_out66.txt` and `c_out67.txt`) and the legacy *FORTRAN* output for this case (the files named, respectively, `fort.66` and `fort.67` in the `test_data/trapping/` folder). Consider that some of the output values will be affected by numeric noise and take substantially different values. However, this is expected for results whose order of magnitude is clearly below the precision level of the calculation, as they represent results appraching zero that were just approximated with different precision.

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
