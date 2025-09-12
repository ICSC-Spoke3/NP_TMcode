# Folder instructions

This directory contains the material to build the project inline documentation with *doxygen*.

## Instructions

The project documentation is managed by *doxygen*, a documentation generator that is able to extract documents directly from properly formatted comment sections of the source code. To build a local instance of project documents, make sure that you have *doxygen* installed, then `cd` into the document source folder (the folder containing the `config.dox` file, specifically `np_tmcode/doc/src`) and finally run:

> doxygen config.dox

*doxygen* will automatically build the HTML structure to cover all the documented source code and it will additionally provide the fundamental structure to prepare a LaTeX formatted version of the documents. These two outputs will be placed, respectively, under the folders `np_tmcode/doc/build/html` and `np_tmcode/doc/build/latex`. If the host system has a recommended LaTeX installation, the full PDF documentation can be built with a further execution of `make` from the `np_tmcode/doc/build/latex` folder.

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
