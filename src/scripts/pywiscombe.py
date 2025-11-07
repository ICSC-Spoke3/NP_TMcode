#!/usr/bin/env python3

#   Copyright (C) 2025   INAF - Osservatorio Astronomico di Cagliari
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

## @package pywiscombe
#  \brief Script to obtain a quick guess of the correct orders for a model.
#
#  This script uses the Wiscombe criterion and its modifications to suggest
#  the proper multipole truncation orders for a particle. It can be used to
#  obtain an estimate of the proper internal order LI (used with the model
#  components or with the model in far field approximation) and the proper
#  external order LE (used for the model in near field approximation).
#
#  The script execution requires python3.

import math
from sys import argv

## \brief Main execution code
#
# `main()` is the function that handles the creation of the script configuration
# and the execution of the comparison. It returns an integer value corresponding
# to the number of detected error-level inconsistencies.
#
# \returns errors: `int` Number of runtime errors (0 means successful run).
def main():
    config = {}
    errors = 0
    try:
        config = parse_arguments()
        if config['help_mode']:
            print_help()
        else:
            if (config['mode'] == ''):
                print("ERROR: no calculation mode chosen (must be either \"--li\" or \"--le\").")
                errors += 1
            if (config['wavelength'] == 0.0):
                print("ERROR: no wavelength specified (missing \"--wave=WAVELENGTH\" option).")
                errors += 1
            if (config['radius'] <= 0.0):
                print("ERROR: invalid particle radius (must be \"--rad=RADIUS\" with RADIUS > 0).")
                errors += 1
            if (errors == 0):
                epsilon = math.sqrt(config['refraction'])
                x = 2 * math.pi * epsilon * config['radius'] / config['wavelength']
                lmax = 0
                if (config['mode'] == 'LI'):
                    lmax = 2 + int(math.ceil(x + 4.05 * math.pow(x, 1.0 / 3.0)))
                elif (config['mode'] == 'LE'):
                    lmax = 1 + int(math.ceil(x + 11.0 * math.pow(x, 1.0 / 3.0)))
                print("Suggested truncation order is Lmax = %d"%lmax)
    except ValueError as ex:
        print(ex)
        print("\nRun \"pywiscombe.py --help\" to get more detailed help.")
        errors = 1
    return errors

## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of mandatory and optional
#  arguments. Mandatory arguments are those required to execute a meaningful
#  parsing and they are limited to the names of the files that need to be read
#  and written. The other arguments affect the format and the data fields that
#  are sought for.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'mode': '',
        'wavelength': 0.0,
        'radius': 0.0,
        'refraction' : 1.0,
        'help_mode': False
    }
    arg_index = 1
    skip_arg = False
    for arg in argv[1:]:
        if skip_arg:
            skip_arg = False
            continue
        split_arg = arg.split('=')
        if (arg.startswith("--li")):
            config['mode'] = 'LI'
        elif (arg.startswith("--le")):
            config['mode'] = 'LE'
        elif (arg.startswith("--wave=")):
            config['wavelength'] = float(split_arg[1])
        elif (arg.startswith("--rad=")):
            config['radius'] = float(split_arg[1])
        elif (arg.startswith("--n=")):
            config['refraction'] = float(split_arg[1])
        elif (arg.startswith("--help")):
            config['help_mode'] = True
        else:
            raise ValueError("Unrecognized argument \'{0:s}\'".format(arg))
        arg_index += 1
    return config

## \brief Print a command-line help summary.
def print_help():
    print("                                                            ")
    print("***                  WISCOMBE CRITERION                  ***")
    print("                                                            ")
    print("Provide a guess of the proper truncation orders for a model.")
    print("                                                                                   ")
    print("Usage: \"./pywiscombe.py --li|le --wave=WAVELENGTH --rad=RADIUS\"                  ")
    print("                                                                                   ")
    print("Valid options are:                                                                 ")
    print("--li                     Run the calculation for an internal field expansion order.")
    print("--le                     Run the calculation for an external field expansion order.")
    print("--wave=WAVELENGTH        Radiation wavelength in meters (mandatory).               ")
    print("--rad=RADIUS             Particle radius in meters (mandatory).                    ")
    print("--help                   Print this help and exit.                                 ")
    print("                                                                                   ")

# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
exit(res)
