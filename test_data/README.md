# Folder instructions

This directory collects the input and output data files to test the code.

## Instructions

The execution of the original code can be controlled through a set of configuration files that define the characteristics of the problem and affect the type of output. In order to allow for consistency tests between the results of the legacy _FORTRAN_ implementation and the new _C++_ version, a set or pre-computed test cases is provided. The test cases are arranged in folders named after the code branch they are intended to test. Namely:

- the `cluster` folder contains input and output data for test cases using clusters of spheres;
- the `sphere` folder contains input and output data for test cases based on a single sphere;
- the `trapping` folder contains input and output data for test cases applying to the `trapping` code.

The legacy _FORTRAN_ code is distributed together with the `np_tmcode` software package, so that the user can verrify the correspondence between results obtained by different code builds. However, some of the _FORTRAN_ based calculations are time consuming. For this reason, the set of pre-computed output files should serve as a first choice for reference.

In general, input data files are identified by names that start with the `D` character, while output files have names that start with the `O` character. The _C++_ implementation of the code prepends a `c_` prefix to all of its output files (while still using the same input ones).

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
