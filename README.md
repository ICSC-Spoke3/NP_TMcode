# NP_TMcode project

This project is the T-matrix based suite of codes by the Messina group (Borghese, Denti, Saija, Iatì) to compute scattering and extinction properties of realistic particles that can be represented as an arbitrary assembly of individually spherically symmetric subparticles, whose classical optical constants are known. Using an arbitrarily large set of arbitrarily small spherically symmetric subparticles one can obviously approximate any real particle, very much as small enough pixels can approximate an image.

The aim of the project, funded by PNRR-CNS, is to refactor the original, very old legacy Fortran codes, casting them in a modern language that allows them to be parallelised and run efficiently on current and forthcoming HPC architectures.

The current implementation offers a set of elementary tests to check that the original FORTRAN code can be compiled and executed on a limited set of pre-defined input data. The functionality of this initial stage can be verified by cloning the gitLab repository on a local machine and building the binaries from the `src` folder.

Distributing the code and its sources is possible under the terms of the GNU GPLv3 license (see *License* below). Use of this code and of any derived implementation should credit the original authors. Scientific publications should do so by citing the following references:

- Saija et al. 2001, ApJ, 559, 993, DOI:10.1086/322350
- Borghese, Denti, Saija 2007, Scattering from Model Nonspherical Particles (ISBN 978-3-540-37413-8), DOI:10.1007/978-3-540-37414-5

*NOTE:* The building process requires a working installation of a C++ and a FORTRAN compiler. Many solutions are available, but the recommended option is the *GNU Compiler Collection* `gcc` with the addition of `g++` and `gfortran`. Single-workstation multi-threaded parallelism is supported through _OpenMP_, while the multi-node code implementation further requires the use of parallel compilers complying with the _MPI_ standard (_OpenMPI_, _MPICH_).

# Acknowledgments

Supported by Italian Research Center on High Performance Computing Big Data and Quantum Computing (ICSC), project funded by European Union - NextGenerationEU - and National Recovery and Resilience Plan (NRRP) - Mission 4 Component 2 within the activities of Spoke 3 (Astrophysics and Cosmos Observations). This work was completed in part at the CINECA GPU HACKATHON 2024, part of the Open Hackathons program. The authors would like to acknowledge OpenACC-Standard.org for their support.

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
