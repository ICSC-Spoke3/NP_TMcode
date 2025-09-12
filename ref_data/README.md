# Folder instructions

This directory contains reference data for model configuration, including, in particular, the dielectric constant files.

## Instructions

Dielectric constant files must be formatted as comma separated values, following the scheme provided by the available example files. `NPtm_code` can extract the necessary dielectric constants from these files, either by interpolating them on an arbitrary set of wavelengths, or just computing the scattering process exactly on the grid of available constants.

If the calculation uses the exact grid matching, then all the dielectric constants *must* be provided on the same set of wavelengths. This constraint does not apply to the interpolation mode, because the dielectric constants will be internally interpolated on the same grid of wavelengths. In this second case, however, it is user's responsibility that the requested wavelength range is covered by all the necessary dielectric constant data.

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
