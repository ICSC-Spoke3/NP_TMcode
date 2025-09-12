# Folder instructions

This directory contains scripts and tools to evaluate the code functionality.

## Instructions

The code migration stage can be considered successfully fulfilled with the solution of the single sphere and of the cluster of spheres cases in *C++*. To test the reliability of the code, the *C++* version needs to produce consistent output with respect to the original *FORTRAN* code. Since the output files are generally too large for manual inspection and they can be affected by different approximations or numeric noise, a series of scripts designed to compare the two versions and to pin-point possible differences is provided in this folder. The scripts are written in *python*, therefore they need an available *python* environment to work.

## Comparing the results of *FORTRAN* and *C++* codes

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

## Estimating the execution time from a terminal log

Performance estimates can be obtained from the code logging system, assuming the uses chose to save the terminal output to some log file. To obtain the time spent by the code in performing a specific operation, the syntax is:

   > $PATH_TO_SCRIPT/pytiming.py --logname=LOG_FILE [--filter=FILTER --threads=NUM_THREADS]

where `LOG_FILE` must be the name of a file containing the output that would normally go to terminal, `FILTER` must be the starting character sequence of the log line identifying the operation that should be taken into account, and `NUM_THREADS` is the number of processes that were used to perform the whole calculation loop. In case no filter is given, the script provides an estimate of the total amount of time spent in doing the calculation. This estimate, however, is known not to be reliable, because it ignores thread concurrency effects. A more accurate estimate of the total time spent in executing the code is always saved in a file named `c_timing.log`.

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
