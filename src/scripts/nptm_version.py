#   Copyright (C) 2026   INAF - Osservatorio Astronomico di Cagliari
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   A copy of the GNU General Public License is distributed along with
#   this program in the COPYING file. If not, see: <https://www.gnu.org/licenses/>.

## @package nptm_version
#  \ingroup python_tools
#  \brief Version of the package, read from `CITATION.cff`.
#
#  The "version" key of `CITATION.cff`, in the root folder of the
#  repository, is the version of the whole package: `configure` passes it to
#  the C++ programs, this module to the Python scripts, and the manual reads
#  it when it is compiled. A script copied out of the repository reports an
#  unknown version.

import os
import re

## \brief Read the version from `CITATION.cff`.
#
#  \param citation_file: `str` Path of the citation file (optional; default:
#  `CITATION.cff` two folders above this module).
#  \return version: `str` The version, or "unknown" if it cannot be read.
def read_version(citation_file=None):
    if (citation_file is None):
        here = os.path.dirname(os.path.abspath(__file__))
        citation_file = os.path.join(here, "..", "..", "CITATION.cff")
    try:
        with open(citation_file, "r") as cff:
            for line in cff:
                match = re.match(r'^version:\s*["\']?([^"\'\s]+)', line)
                if (match is not None):
                    return match.group(1)
    except OSError:
        pass
    return "unknown"

## \brief Version of the package.
__version__ = read_version()
