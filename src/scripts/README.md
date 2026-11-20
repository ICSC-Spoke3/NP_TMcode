# Folder instructions

This directory contains scripts and tools to evaluate the code functionality.

## Instructions

The code migration stage can be considered successfully fulfilled with the solution of the single sphere and of the cluster of spheres cases in *C++*. To test the reliability of the code, the *C++* version needs to produce consistent output with respect to the original *FORTRAN* code. Since the output files are generally too large for manual inspection and they can be affected by different approximations or numeric noise, a series of scripts designed to compare the two versions and to pin-point possible differences is provided in this folder. The scripts are written in *python*, therefore they need an available *python* environment to work.

The scripts print their version with `--version`: the `version` key of `CITATION.cff`, read by the `nptm_version.py` module (`unknown` for a script copied out of the repository).

## Create models from YAML input

The execution of `NP_TMcode` applications requires machine readable input files to describe the particle model and the radiation field properties. This files are ASCII files that can be edited with common text editors, but they obey a strict format which complies with the requirements of the original code implementation. In order to make the model creation easier, `NP_TMcode` provides a model editor script named `model_maker.py`. This script reads a model description from a `YAML` file and creates the necessary input configuration for the `NP_TMcode` applications. To run the script, just issue:

   > $PATH_TO_SCRIPT/model_maker.py YAML_CONFIGURATION_FILE

where `YAML_CONFIGURATION_FILE` must be a valid `YAML` model description (examples are given in the `test_data` folder). More details are given issuing:

   > $PATH_TO_SCRIPT/model_maker.py --help

## Compare the results of *FORTRAN* and *C++* codes

1. Follow the instructions to build and run the *FORTRAN* and *C++* versions of the code.
2. Run the `pycompare.py` script providing the following minimal arguments:

   > $PATH_TO_SCRIPT/pycompare.py --ffile FORTRAN_OUTPUT --cfile C++_OUTPUT

(The above assumes that `PATH_TO_SCRIPT` is a variable that expands to the path were the script is located). The required output files are called by default `OSPH` and `OCLU` by *FORTRAN* and `c_OSPH` and `c_OCLU` by *C++*, depending on whether the sphere or the cluster case was executed. 

3. Check the output of the script to verify that it detects 0 errors and finishes in a `SUCCESS` state.
4. In case of necessity, add the `--html` argument to produce an *HTML* log showing the possible differences and a classification of their severity.
5. Issuing:

   > $PATH_TO_SCRIPT/pycompare.py --help

or just:

   > $PATH_TO_SCRIPT/pycompare.py

will print a help screen, giving a brief explanation of all the possible options.

## Estimate the recommended orders for running a model

This script applies the Wiscombe criterion and its modifications to estimate the recommended internal and external orders to run a model. To use it, issue:

   > $PATH_TO_SCRIPT/pywiscombe.py --li|le --wave=WAVELENGTH --rad=RADIUS

Wavelength and radius must be expressed in meters, with the radius being equal to the minimum radius encircling the particle. More details are given with:

   > $PATH_TO_SCRIPT/pywiscombe.py --help

## Estimate the memory needed by a model

`guess_resources.py` reads the input files of a model and prints the particle radius and the memory estimates that the programs log and check before a calculation (host memory per thread and per MPI process, output data, GPU memory of the matrix inversion). The runtime options of the geometry file (`USE_REFINEMENT`, `INVERSION`, `DEVICE_POLICY`, `HOST_RAM_GB`, `GPU_RAM_GB`) are taken into account; the build and the run are described by options:

   > $PATH_TO_SCRIPT/guess_resources.py --edfb EDFB_FILE --geom GEOM_FILE [--backend magma|cublas|lapack|lucin] [--threads N] [--processes N] [--gpus N] [--host-ram GIB] [--gpu-ram GIB]

The estimates are computed by the `nptm_estimates.py` module, which `model_maker.py` also uses. They reproduce the `get_size()` functions of the C++ code: a change in those functions must be copied to the module. The "minimum encircling radius" printed by both scripts is the one used by the code to choose the external order: the radius of the smallest sphere centred in the origin of the coordinates that contains the particle (the external sphere, for inclusions). This is the radius to give to `pywiscombe.py --le`.

## Compute the mass, the tensor of inertia and the filling factor of a model

`inertia.py` (experimental) reads the scatterer and geometry files of a model and computes the mass of each sphere type and of the particle, the tensor of inertia with respect to the origin and to the centre of mass, and the filling factor: the volume of the particle over the volume of the homogeneous ellipsoid with the same mass and principal moments of inertia.

   > $PATH_TO_SCRIPT/inertia.py --scat EDFB_FILE --geom GEOM_FILE [--spew=W1,W2,...] [--medium=W]

The materials are the distinct dielectric constants of the first wavelength, in the order of the scatterer file (the script lists them); `--spew` gives their specific weights in g/cm<sup>3</sup> (one value for all of them, default 1). Transition layers change density with the profile that the code uses for the dielectric constant; `--medium` is the specific weight of the external medium, which a last transition layer reaches (default 0). For `np_inclusion` models, the external sphere is centred in the origin and filled with its material around the spheres.

## Estimate the execution time from a terminal log

Performance estimates can be obtained from the code logging system, assuming the uses chose to save the terminal output to some log file. To obtain the time spent by the code in performing a specific operation, the syntax is:

   > $PATH_TO_SCRIPT/pytiming.py --logname=LOG_FILE [--filter=FILTER --threads=NUM_THREADS]

where `LOG_FILE` must be the name of a file containing the output that would normally go to terminal, `FILTER` must be the starting character sequence of the log line identifying the operation that should be taken into account, and `NUM_THREADS` is the number of processes that were used to perform the whole calculation loop. In case no filter is given, the script provides an estimate of the total amount of time spent in doing the calculation. This estimate, however, is known not to be reliable, because it ignores thread concurrency effects. A more accurate estimate of the total time spent in executing the code is always saved in a file named `c_timing.log`.

## Extract parameters from the code output

The legacy and `HDF5` outputs of `NP_TMcode` applications contain all the values computed by the code. However, in many cases the user is interested in extracting plots of cross-sections, forces and torques as functions of wavelength. The `parse_output.py` script performs this task, on either output: the legacy file (`c_OSPH`, `c_OCLU`, `c_OINCLU`, or the output of the FORTRAN codes) or the `HDF5` one (`c_OSPH.hd5`, `c_OCLU.hd5`, `c_OINCLU.hd5`), which give the same `CSV` files. The format is recognized from the file; reading `HDF5` files requires the `h5py` Python module. To use the script after running a model, issue:

   > $PATH_TO_SCRIPT/parse_output.py --in CODE_RESULT_FILE --out OUTPUT_BASE_NAME --app=[SPH|CLU|INCLU]

This will parse the output of `np_sphere` (`app=SPH`), `np_cluster` (`app=CLU`) or `np_inclusion` (`app=INCLU`) and put the results in `CSV` text files, named after the given `OUTPUT_BASE_NAME`, but with distinguishing suffixes such as `_ics` for integrated cross-sections, `_dcs` for differential ones, `_irp` for integrated radiation pressure forces, etc. The files of `np_cluster` and `np_inclusion` come in pairs, for the polarization states -1 (`_ics1`, ...) and +1 (`_ics2`, ...); the `_drp` files have a row for each direction in which the code computes forces and torques (every direction with `ISAM` < 0, the first scattered direction of each incident direction otherwise). The wavelength is 2 * PI / VK, or 2 * PI / (XI * VK) when XI is a scale factor for lengths (`IDFC` < 0). More details are given issuing:

   > $PATH_TO_SCRIPT/parse_output.py --help

## Extract the dynamic range of a complex matrix

This script examines a complex matrix and extracts histograms of the orders of magnitude of the real parts, the imaginary parts and the absolute values. Optionally, if `matplotlib` is available, a plot of the dynamic range is drawn. To use this script, issue:

   > $PATH_TO_SCRIPT/pydynrange.py MATRIX_FILE

where `MATRIX_FILE` must be a `CSV` representation of the matrix. More details are given issuing:

   > $PATH_TO_SCRIPT/pydynrange.py --help

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
