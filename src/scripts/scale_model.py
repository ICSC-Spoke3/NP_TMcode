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

## @package scale_model
#  \brief Script to re-scale an existing model by multiplicative factor.
#
#  This script is designed to create scaled versions of pre-existing
#  models, to allow for quick production of sets of particle models
#  characterized by the same structure, but with a different size.
#
#  The script execution requires python3.

# import pdb # uncomment this line for debugging
import re
from sys import argv

## \cond
int_reg = re.compile(r'-?[0-9]+')
float_reg = re.compile(r'-?[0-9]+\.[0-9]+([eEdD][-+]?[0-9]+)?')
## \endcond

## \brief Main execution code
#
# `main()` is the function that handles the creation of the script configuration
# and the execution of the comparison. It returns an integer value corresponding
# to the script exit code.
#
# \returns errors: `int` Number of runtime errors (0 means successful run).
def main():
    config = {}
    errors = 0
    try:
        config = parse_arguments()
    except ValueError as ex:
        print(ex)
        print("\nType \"scale_model.py --help\" to get more detailed help.")
        errors = 1
    if config['help_mode']:
        config['help_mode'] = True
        print_help()
    else:
        if (config['input_name'] == ''):
            print("ERROR: no input file name given (missing \"--in FILE\" option).")
            errors += 1
        if (config['mode'] == ''):
            print("ERROR: no mode specifier given (missing \"--mode=MODE\" option).")
            errors += 1
        if (errors == 0):
            if (config['mode'] == 'edfb'):
                try:
                    errors += scale_legacy_edfb(config)
                except Exception as ex:
                    print(ex)
                    errors += 1
            elif (config['mode'] == 'geom'):
                try:
                    errors += scale_legacy_geom(config)
                except Exception as ex:
                    print(ex)
                    errors += 1
            else:
                print("ERROR: unknown mode \"%s\" (options are edfb|geom)"%(config['mode']))
                errors += 1
    return errors

## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of mandatory and optional
#  arguments. Mandatory arguments are those required to execute a meaningful
#  parsing and they are limited to the names of the files that need to be read
#  and the operational mode ("edfb" for scattering model file, "geom" for
#  geometry model file).
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'input_name': '',
        'output_name': '',
        'mode': '',
        'scale': 1.0,
        'help_mode': False
    }
    arg_index = 1
    skip_arg = False
    for arg in argv[1:]:
        if skip_arg:
            skip_arg = False
            continue
        split_arg = arg.split('=')
        if (arg.startswith("--in")):
            config['input_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--out")):
            config['output_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--mode=")):
            config['mode'] = split_arg[1]
        elif (arg.startswith("--scale=")):
            config['scale'] = float(split_arg[1])
        elif (arg.startswith("--help")):
            config['help_mode'] = True
        else:
            raise ValueError("Unrecognized argument \'{0:s}\'".format(arg))
        arg_index += 1
    if (config['output_name'] == ''):
        config['output_name'] = config['input_name'] + '_scaled'
    return config

## \brief Scale a model legacy ScattererConfiguration file.
#
#  \param config: `dict` A dictionary containing the script configuration.
#  \return errors: `int` The number of encountered errors.
def scale_legacy_edfb(config):
    errors = 0
    input_file = open(config['input_name'], 'r')
    output_file = open(config['output_name'], 'w')
    file_line = input_file.readline() # NSPH IES
    file_line = file_line.replace("D", "E").replace("d", "e")
    configurations = 0
    iter_numbers = int_reg.finditer(file_line)
    n_groups = []
    for ni in iter_numbers:
        n_groups.append(ni.group())
    nsph = int(n_groups[0])
    ies = int(n_groups[1])
    output_file.write(file_line)
    file_line = input_file.readline() # EXDC WP XIP IDFC NXI INSTPC INSN
    file_line = file_line.replace("D", "E").replace("d", "e")
    configurations = 0
    iter_numbers = float_reg.finditer(file_line)
    n_ends = []
    for ni in iter_numbers:
        n_ends.append(ni.end())
    iter_numbers = int_reg.finditer(file_line[n_ends[2]:])
    n_groups = []
    for ni in iter_numbers:
        n_groups.append(ni.group())
    nxi = int(n_groups[1])
    instpc = int(n_groups[2])
    output_file.write(file_line)
    # Read and copy the wavelength scales
    if (instpc == 1):
        file_line = input_file.readline()
        file_line = file_line.replace("D", "E").replace("d", "e")
        output_file.write(file_line)
    else:
        for fli in range(nxi):
            file_line = input_file.readline()
            file_line = file_line.replace("D", "E").replace("d", "e")
            output_file.write(file_line)
    # Read and copy the vector of sphere IDs
    vec_ids = []
    while (len(vec_ids) < nsph):
        file_line = input_file.readline()
        file_line = file_line.replace("D", "E").replace("d", "e")
        iter_numbers = int_reg.finditer(file_line)
        for ni in iter_numbers:
            vec_ids.append(int(ni.group()))
        output_file.write(file_line)
    configurations = 0
    for ci in range(len(vec_ids)):
        if (vec_ids[ci] > ci):
            configurations += 1
    # Read and scale the configuration size
    for ci in range(configurations):
        # work on a configuration ...
        file_line = input_file.readline()
        file_line = file_line.replace("D", "E").replace("d", "e")
        iter_numbers = int_reg.finditer(file_line)
        n_groups = []
        for ni in iter_numbers:
            n_groups.append(ni.group())
        num_layers = int(n_groups[0])
        iter_numbers = float_reg.finditer(file_line)
        n_groups = []
        for ni in iter_numbers:
            n_groups.append(ni.group())
        r_sph = float(n_groups[0]) * config['scale']
        str_line = "{0:3d}   {1:14.6e}\n".format(num_layers, r_sph)
        output_file.write(str_line)
        for cj in range(num_layers):
            file_line = input_file.readline()
            file_line = file_line.replace("D", "E").replace("d", "e")
            output_file.write(file_line)
    # Read and copy the optical constants
    file_line = input_file.readline()
    while(file_line != ""):
        file_line = file_line.replace("D", "E").replace("d", "e")
        output_file.write(file_line)
        file_line = input_file.readline()
    input_file.close()
    output_file.close()
    return errors

## \brief Scale a model legacy GeometryConfiguration file.
#
#  \param config: `dict` A dictionary containing the script configuration.
#  \return errors: `int` The number of encountered errors.
def scale_legacy_geom(config):
    errors = 0
    input_file = open(config['input_name'], 'r')
    output_file = open(config['output_name'], 'w')
    file_line = input_file.readline() # HEADER LINE
    iter_numbers = int_reg.finditer(file_line)
    n_groups = []
    for ni in iter_numbers:
        n_groups.append(ni.group())
    nsph = int(n_groups[0])
    output_file.write(file_line)
    for si in range(nsph):
        file_line = input_file.readline() # First data line
        file_line = file_line.replace("D", "E").replace("d", "e")
        iter_numbers = float_reg.finditer(file_line)
        n_groups = []
        for ni in iter_numbers:
            n_groups.append(ni.group())
        if (len(n_groups) == 3):
            # do it
            sph_x = float(n_groups[0]) * config['scale']
            sph_y = float(n_groups[1]) * config['scale']
            sph_z = float(n_groups[2]) * config['scale']
            str_line = "   {0:15.7e}   {1:15.7e}   {2:15.7e}\n".format(
                sph_x, sph_y, sph_z
            )
            output_file.write(str_line)
        else:
            print("ERROR: sphere coordinates vectors not in place!")
            errors += 1
            break # si loop on spheres
    # Read and parse the directional settings
    for li in range(2):
        file_line = input_file.readline()
        file_line = file_line.replace("D", "E").replace("d", "e")
        output_file.write(file_line)
    # Read and copy runtime options
    file_line = input_file.readline()
    while (file_line != ""):
        output_file.write(file_line)
        file_line = input_file.readline()
    input_file.close()
    output_file.close()
    return errors

## \brief Print a command-line help summary.
def print_help():
    print("                                                  ")
    print("***                SCALE MODEL                 ***")
    print("                                                  ")
    print("Scale model input files by multiplicative factor. ")
    print("                                                                           ")
    print("Usage: \"./scale_model.py --in INPUT --mode=MODE [--out OUTPUT] [OPTIONS]\"")
    print("                                                                           ")
    print("Valid options are:                                                         ")
    print("--in INPUT               Original model to be scaled (mandatory).          ")
    print("--out OUTPUT             Name for the rescaled model (optional).           ")
    print("--mode=[edfb|geom]       Type of input to be processed (mandatory).        ")
    print("--scale=SCALE            Scale to be applied (optional, default is 1).     ")
    print("--help                   Print this help and exit.                         ")
    print("                                                                           ")

# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
exit(res)
