# Folder instructions

This directory contains the results of example calculations executed on particle models characterized by high complexity in terms of number of particle constituents, materials and spherical monomer types.

## Instructions

This folder contains the input parameters and the code output for a set of particle models, illustrated in the 3D rendering images provided for each structure. Each of the proposed structures is solved for different materials, after scaling the particle structure by various factors, thus effectively representing particles of different size. The results of the calculations are stored in sub-folders divided by the material of the particle and, then, by the overall size of the aggregate.

## Folder contents

Each model is characterized by a root file name encoded as follows:

nXXX_stX_rmXXum_rndXX_[scX_Y]_DESC

where the variable parts of the file names indicate, respectively:

- nXXX: the total number of spheres used by the model
- stX: the number of sphere types defined in the model
- rmXXum: the particle maximum radius in micrometers at scale 1
- scX_Y: if present, indicates that the model is rescaled by factor X.Y
- DESC: a file content description which can take the values
  - config.yml: the model configuration YAML file
  - model.png: a rendering of the particle structure
  - .txt: particle description parameters (metadata)
  - DCLU: model geometry configuration file (NP_TMcode input)
  - DEDFB: model scatterer configuration file (NP_TMcode input)
  - irp.csv: comma-separated value data for integrated radiation pressure
  - ics.csv: comma-separated value data for integrated cross-sections

# License

   Copyright (C) 2026   INAF - Osservatorio Astronomico di Cagliari

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
