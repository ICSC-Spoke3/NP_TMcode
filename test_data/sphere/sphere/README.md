# Folder instructions

This directory contains test data for the single sphere case.

## Instructions

`sph` and its _C++_ implementation, `np_sphere`, are designed to perform the simplest case calculation, namely the scattering of incident radiation on a single sphere. To perform the calculation, the two following formatted files need to be provided:

- DEDFB
```
   SPHERE_NUMBER   EXTERNAL_SPHERE_FLAG
   EXT_DIEL_CONST   PEAK_WAVE   PEAK_SCALE   DIEL_TYPE   SCALE_NUMBER   SCALE_STEP_FLAG   CONTROL_VARIABLE
   SCALE_VECTOR   [SCALE_VECTOR_STEP]
        |
	v
   SPHERE_ID_VECTOR ->
   LAYER_NUMBER_SPHERE_1   SPHERE_1_RADIUS
   FRACTIONAL_BREAK_RADII
             |
	     v
   DIEL_CONST_VECTOR
          |
	  v
   EOF_CODE
```
where the different lines have the following roles:
1. declare the number of spheres and whether to add an external layer
2. define the external dielectric constant and the scaling configuration
3. define the vector of scales (either explicitly, with one element per row, or in steps, with only first element and step declared)
4. create a vector of IDs to attach to the spheres
5. define the vectors specifying the number of layers and the radii of the spheres (one sphere per line)
6. define the fractional break radii (one line per layer)
7. define the vector of dielectric constants or starting functions (one line per layer in sphere)
8. an end-of-file code (generally 0)

- DSPH
```
  SPHERE_NUMBER  MAXIMUM_L_ORDER  POLARIZATION_STATUS  TRANSITION_SHARPNESS_1  TRANSITION_SHARPNESS_2  GEOMETRY
  STARTING_INC_THETA  INC_THETA_STEP  FINAL_INC_THETA  STARTING_SCA_THETA  SCA_THETA_STEP  FINAL_SCA_THETA
  STARTING_INC_PHI  INC_PHI_STEP  FINAL_INC_PHI  STARTING_SCA_PHI  SCA_PHI_STEP  FINAL_SCA_PHI
  OUTPUT_T-MATRIX_SCALE_NUMBER
EOF_CODE
```
where the different lines have the following roles:
1. general configuration of the scattering problem, with some specification of the transition between materials
2. definition of the elevation angle arrays for the incident and scattered radiation fields
3. definition of the azimuth angle arrays for the incident and scattered radiation fields
4. the number of the calculation scale for which the T-matrix will be saved to file
5. an end-of-file code (generally 0)

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
