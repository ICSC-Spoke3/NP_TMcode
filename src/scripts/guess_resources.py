#!/usr/bin/env python3

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

## @package guess_resources
#  \brief Script to estimate the computational cost of a model.
#
#  Script to assist in the creation of model input files starting
#  from a YAML descriptor.
#
#  The script requires python3.

import math
import numpy as np
import pdb
import re
from sys import argv

## \cond
__version__ = "0.10.10"
int_reg = re.compile(r'[-+]?[0-9]+')
number_reg = re.compile(r'[-+]?[0-9]+\.[0-9]+([eE][-+]?[0-9]+)?')
## \endcond

## \brief Main execution code.
#
# `main()` is the function that handles the creation of the code configuration.
# It returns an integer value as exit code, using 0 to signal successful execution.
#
# \returns result: `int` Exit code (0 = SUCCESS).
def main():
    result = 0
    config = parse_arguments()
    if config['version_mode']:
        print("guess_resources.py v%s."%__version__)
    elif (config['help_mode']):
        print_help()
    else:
        if (config['edfb_file'] == ''):
            print("ERROR: missing scatterer configuration file (--edfb EDFB_FILE).")
            result = 1
        elif (config['geom_file'] == ''):
            print("ERROR: missing geometry configuration file (--geom GEOM_FILE).")
            result = 1
        else:
            model = scan_model(config['edfb_file'], config['geom_file'])
            if (model is not None):
                print_model_summary(model)
                print_resource_summary(model)
            else:
                result = 1
    return result

## \brief Transform a callable iterator object into an array.
#
#  \param[in,out] my_iter: `callable_iterator` Iterator to be transformed (gets consumed).
#  \returns array: `array`-like An array of items from the iterator.
def iter_to_array(my_iter):
    array = []
    for item in my_iter:
        array.append(item)
    return array

## \brief Compute the particle's filling factor.
#
#  In this application the filling factor is defined as the ratio of the
#  sum of the volumes of the spheres that compose the particle with respect
#  to the volume of a homogeneous ellipsoid that has the same rotational
#  properties as the particle. The result provided by this script is an
#  approximation that assumes the same specific weight for all monomers.
#  More accurate estimates (including the actual mass of the particle)
#  should be computed via the `inertia.py` script.
#
#  \param model: `dict` Model description dictionary.
#  \return ff: `float` The particle's filling factor.
def filling_factor(model):
    ff = 0.0
    N = model['nsph']
    Itot = np.zeros((3, 3))
    #vtot = get_total_volume(stypes, config)
    vtot = 0.0
    volumes = np.zeros((N))
    gf = 4.0 * math.pi / 3.0
    for i in range(model['nsph']):
        sph_type_index = model['vec_types'][i] - 1
        ros = model['ros'][sph_type_index]
        ros3 = ros * ros * ros
        volumes[i] = gf * ros3
        vtot += volumes[i]
    # for i ends here
    masses = volumes * 1000.0
    mtot = np.sum(masses, axis=0)

    positions = np.array([
        [model['vec_x'][i], model['vec_y'][i], model['vec_z'][i]] for i in range(N)
    ])
    center_of_mass = np.sum(masses[:, np.newaxis] * positions, axis=0) / mtot
    for i in range(N):
        M = masses[i]
        typeID = model['vec_types'][i] - 1
        R = model['ros'][typeID]
        x = model['vec_x'][i] - center_of_mass[0]
        y = model['vec_y'][i] - center_of_mass[1]
        z = model['vec_z'][i] - center_of_mass[2]
        
        # 1. Get the inertia tensor of each sphere
        I_sph_val = (2/5) * M * R**2
        I_sph = np.eye(3) * I_sph_val
        
        # 2. Get the transfer tensor through the parallel axes theorem.
        # Diagonals: M * (sum of squared other coordinates)
        # Off-diagonals: -M * (product of coordinates)
        I_transfer = M * np.array([
            [y**2 + z**2, -x * y, -x * z],
            [-y * x, x**2 + z**2, -y * z],
            [-z * x, -z * y, x**2 + y**2]
        ])
        Itot += (I_sph + I_transfer)
    # i loop ends here
    
    eigenval, eigenvec = np.linalg.eigh(Itot)
    print("INFO: ellipsoid directions:")
    for vi in range(3): print(" "*5, eigenvec[vi])
    gf = 5.0 / (2.0 * mtot) 
    I1, I2, I3 = eigenval
    a = np.sqrt(max(0.0, gf * (I2 + I3 - I1)))
    b = np.sqrt(max(0.0, gf * (I1 + I3 - I2)))
    c = np.sqrt(max(0.0, gf * (I1 + I2 - I3)))
    # Volume of the ellipsoid that has the same moments of inertia as the particle
    vell = 4.0 * math.pi * a * b * c / 3.0
    ff = vtot / vell
    return ff

## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of optional arguments.
#  The purpose of this function is to parse the command line in search for
#  such arguments and prepare the execution accordingly.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'edfb_file': '',
        'geom_file': '',
        'help_mode': False,
        'version_mode': False
    }
    skip_arg = False
    dict_key = ''
    for arg in argv[1:]:
        if (skip_arg):
            if (dict_key != ''):
                config[dict_key] = arg
                dict_key = ''
                skip_arg = False
                continue
        if (arg.startswith("--help")):
            config['help_mode'] = True
        elif (arg.startswith("--version")):
            config['version_mode'] = True
        elif (arg.startswith("--edfb")):
            dict_key = 'edfb_file'
            skip_arg = True
        elif (arg.startswith("--geom")):
            dict_key = 'geom_file'
            skip_arg = True
        else:
            raise Exception("Unrecognized argument \'{0:s}\'".format(arg))
    return config

## \brief Print a command-line help summary.
def print_help():
    print("###############################################           ")
    print("#                                             #           ")
    print("#          NPTM_code GUESS_RESOURCES          #           ")
    print("#                                             #           ")
    print("###############################################           ")
    print("                                                          ")
    print("Evaluate the resources required by a model.               ")
    print("                                                          ")
    print("Usage: \"./guess_resources.py --edfb EDFB --geom GEOM\"   ")
    print("                                                          ")
    print("EDFB and GEOM must be valid NPTM_code configuration files.")
    print("                                                          ")
    print("Valid options are:                                        ")
    print("--help                Print this help and exit.           ")
    print("--version             Print script version and exit.      ")
    print("                                                          ")

## \brief Print a summary of model properties.
#
#  This function provides a summary of useful information concerning the
#  radii of the particle monomers and of the equivalent mass sphere, to
#  assist in the selection of the proper starting orders.
#
#  \param[in] model: `dict` A model description dictionary.
def print_model_summary(model):
    avgX = 0.0
    avgY = 0.0
    avgZ = 0.0
    avgR = 0.0
    cmX = 0.0
    cmY = 0.0
    cmZ = 0.0
    Rmin = 0.0
    Rmax = 0.0
    Reqm = 0.0
    R3tot = 0.0
    Rcirc = 0.0
    square_farthest = 0.0
    # breakpoint()
    for i in range(model['nsph']):
        sph_type_index = model['vec_types'][i] - 1
        ros = model['ros'][sph_type_index]
        ros3 = ros * ros * ros
        avgX += model['vec_x'][i]
        avgY += model['vec_y'][i]
        avgZ += model['vec_z'][i]
        avgR += ros
        cmX += (ros3 * model['vec_x'][i])
        cmY += (ros3 * model['vec_y'][i])
        cmZ += (ros3 * model['vec_z'][i])
        R3tot += ros3
        if (ros > Rmax):
            Rmax = ros
        if (Rmin == 0.0 or ros < Rmin):
            Rmin = ros
    Reqm = math.pow(R3tot, 1.0 / 3.0)
    avgX /= model['nsph']
    avgY /= model['nsph']
    avgZ /= model['nsph']
    avgR /= model['nsph']
    avgR3 = avgR * avgR * avgR
    cmX /= (avgR3 * model['nsph'])
    cmY /= (avgR3 * model['nsph'])
    cmZ /= (avgR3 * model['nsph'])
    for i in range(model['nsph']):
        sph_type_index = model['vec_types'][i] - 1
        ros = model['ros'][sph_type_index]
        dX = model['vec_x'][i] - avgX
        dY = model['vec_y'][i] - avgY
        dZ = model['vec_z'][i] - avgZ
        square_range = dX * dX + dY * dY + dZ * dZ + ros * ros
        if (square_range > square_farthest):
            square_farthest = square_range
    Rcirc = math.sqrt(square_farthest)
    if (model['app'] == "SPHERE"):
        print("INFO: maximum expansion order LM = %d."%model['li'])
    else:
        print("INFO: maximum internal expansion order LI = %d."%model['li'])
        print("INFO: maximum external expansion order LE = %d."%model['le'])
    if (Rmax == Rmin):
        print("INFO: monomer radius is Rsph = %.5em"%Rmin)
    else:
        print("INFO: smallest monomer radius is Rmin = %.5em"%Rmin)
        print("INFO: largest monomer radius is Rmax = %.5em"%Rmax)
    print("INFO: equivalent volume radius is Reqv = %.5em"%Reqm)
    print("INFO: minimum encircling radius is Rcirc = %.5em"%Rcirc)
    print(
        "INFO: geometric center at [{0:.5e}, {1:.5e}, {2:.5e}]".format(
            avgX, avgY, avgZ
        )
    )
    print(
        "INFO: center of filled volume at [{0:.5e}, {1:.5e}, {2:.5e}]".format(
            cmX, cmY, cmZ
        )
    )
    fill_factor = filling_factor(model)
    print("INFO: the particle's filling factor is %.5g"%fill_factor)

## \brief Print a summary of resources needed per iteration.
#
#  \param[in] model: `dict` A model description dictionary.
def print_resource_summary(model):
    host_mem_gb = 0.0
    gpu_mem_gb = 0.0
    half_layers = 1 + int(model['max_layers'] / 2)
    nhspo = max(model['npnt'], model['npntts']) * 2 - 1
    if (model['app'] == 'SPHERE'):
        nlmmt = 2 * model['li'] * (model['li'] + 2)
        host_mem_gb += 16 * 2 * model['li'] # RMI + REI
        host_mem_gb += 16 * 4 * nlmmt # W
        host_mem_gb += 16 * 21 # FSAS + SAS + VINTS
        host_mem_gb += 8 * 7 # SSCS + SEXS + SABS + SQSCS + SQEXS + SQABS + GCSV
        host_mem_gb += 8 * (1 + model['max_layers']) # ROS + RC
        host_mem_gb += 8 # IOG + NSHL
        host_mem_gb += 16 * 2 * nhspo # RIS + DLRI
        host_mem_gb += 16 * half_layers # DC0
        host_mem_gb += 16 + 8 # VKT + VSZ
        host_mem_gb += 8 * model['max_layers'] # RCF
        host_mem_gb += 8 * (16 + 16) # CMULLR + CMUL
        host_mem_gb += 8 * 3 * 11 # UNIT VECTORS
        host_mem_gb += 8 + 8 # ARGI + ARGS
        host_mem_gb += 16 * 16 # VINT
        host_mem_gb += 8 * (model['li'] * 12 + 2) # ZPV + GAPS
        host_mem_gb += 8 * (2 * 4 + 4 * 4) # TQSE + TQSPE + TQSS + TQSPS
        host_mem_gb += 16 * half_layers # DC0M
        host_mem_gb += 8 * model['nxi'] # XIV
        host_mem_gb += 16 * 3 # ARG + S0 + TFSAS
    elif (model['app'] == 'CLUSTER'):
        nsph = model['nsph']
        ncou = nsph * (nsph - 1)
        litpo = 2 * model['li'] + 1
        litpos = litpo * litpo
        lmtpo = 2 * model['li'] * model['le'] + 1
        lmtpos = lmtpo * lmtpo
        lm = max(model['li'], model['le'])
        lmpo = lm + 1
        ndi = nsph * model['li'] * (model['li'] + 2)
        ndit = ndi + ndi
        nlem = model['le'] * (model['le'] + 2)
        nlemt = nlem + nlem
        nv3j = lm * lmpo * (2 * lm + 7) / 6
        host_mem_gb += 16 * ncou * litpo # VH
        host_mem_gb += 16 * nsph * lmtpo # VJ0
        host_mem_gb += 16 * ncou * litpos # VYHJ
        host_mem_gb += 16 * nsph * lmtpos # VYJ0
        host_mem_gb += 16 * 2 * nsph * model['li'] # RMI + REI
        host_mem_gb += 16 * nlemt * 4 # W
        host_mem_gb += 16 * nlemt * nlemt # AM0M
        host_mem_gb += 16 * nsph * 5 # FSAS + SAS
        host_mem_gb += 16 * nsph * 16 # VINTS
        host_mem_gb += 8 * nv3j # V3J0
        host_mem_gb += 8 * nsph * 11 # sphere CSs and Qs + coords and radii
        host_mem_gb += 8 * nsph * model['max_layers'] # RC
        host_mem_gb += 4 * lmpo * lm # IND3J
        host_mem_gb += 4 * 2 * nsph # IOG + NSHL
        host_mem_gb += 16 * (16 + 16 + 8) # VINT + VINTM + SCSCP + ECSCP + SCSCPM + ECSCPM
        host_mem_gb += 16 * nhspo * 2 # RIS + DLRI
        host_mem_gb += 16 * half_layers # DC0
        host_mem_gb += (16 + 8) * nsph # VKT + VSZ
        host_mem_gb += 16 * 5 # TFSAS + TSAS
        host_mem_gb += 8 * 4 # GCS + SCS + ECS + ACS
        host_mem_gb += 8 * lmtpo # RAC3J
        host_mem_gb += 16 * 2 * ndi * nlem # GIS + GLS
        host_mem_gb += 16 * ndit * nlemt # SAM
        host_mem_gb += 16 * ndit * ndit # AM
        gpu_mem_gb += 16 * ndit * ndit # AM
        host_mem_gb += 8 * nsph * model['max_layers'] # RCF
        host_mem_gb += 8 * (16 + 16 + 16 + 16) # CMULLR + CMUL + CEXTLR + CEXT
        host_mem_gb += 8 * 3 * 11 # UNIT VECTORS
        host_mem_gb += 8 + 8 # ARGI + ARGS
        host_mem_gb += 8 * (lm * 12 + nsph) # ZPV + GAPS
        host_mem_gb += 8 * (3 + 6 + 6) # GAP + GAPV + GAPM
        host_mem_gb += 16 * (6 + 6) # GAPP + GAPPM
        host_mem_gb += 8 * (2 * nsph + 4 * nsph) # TQSE + TQSPE + TQSS + TQSPS
        host_mem_gb += 8 * (2 * 4 + 4 * 4) # TQCE + TQCPE + TQCS + TQCPS
        host_mem_gb += 8 * (3 + 3) # TQEV + TQSV
        host_mem_gb += 16 * nsph * half_layers # DC0M
        host_mem_gb += 8 * model['nxi'] # XIV
        host_mem_gb += 16 * 4 # ARG + S0 + S0M + CCSAM
    elif (model['app'] == "INCLUSION"):
        nsph = model['nsph']
        ncou = nsph * (nsph - 1)
        litpo = 2 * model['li'] + 1
        litpos = litpo * litpo
        lmtpo = 2 * model['li'] * model['le'] + 1
        lmtpos = lmtpo * lmtpo
        lm = max(model['li'], model['le'])
        lmpo = lm + 1
        ndi = nsph * model['li'] * (model['li'] + 2)
        ndit = ndi + ndi
        nlem = model['le'] * (model['le'] + 2)
        nlemt = nlem + nlem
        nv3j = lm * lmpo * (2 * lm + 7) / 6
        host_mem_gb += 16 * ncou * litpo # VH
        host_mem_gb += 16 * nsph * lmtpo # VJ0
        host_mem_gb += 16 * ncou * litpos # VYHJ
        host_mem_gb += 16 * nsph * lmtpos # VYJ0
        host_mem_gb += 16 * 2 * nsph * model['li'] # RMI + REI
        host_mem_gb += 16 * nlemt * 4 # W
        host_mem_gb += 16 * (8 * model['le']) # RM0 + RE0 + RMW + REW + TM + TE + TM0 + TE0
        host_mem_gb += 16 * nlemt * nlemt # AM0M
        host_mem_gb += 8 * nv3j # V3J0
        host_mem_gb += 16 * 3 * 4 # FSAC + SAC + FSACM
        host_mem_gb += 16 * (16 + 16 + 8) # VINT + VINTM + SCSCP + ECSCP + SCSCPM + ECSCPM
        host_mem_gb += 8 * nsph * 4 # sphere coords and radii
        host_mem_gb += 8 * nsph * model['max_layers'] # RC
        host_mem_gb += 4 * lmpo * lm # IND3J
        host_mem_gb += 4 * 2 * nsph # IOG + NSHL
        host_mem_gb += 16 * nhspo * 2 # RIS + DLRI
        host_mem_gb += 16 * (half_layers + 1) # DC0
        host_mem_gb += (16 + 8) * nsph # VKT + VSZ
        host_mem_gb += 8 * lmtpo # RAC3J
        host_mem_gb += 16 * nlemt * ndi # AT
        host_mem_gb += 16 * ndit * ndit # AM
        gpu_mem_gb += 16 * ndit * ndit # AM
        host_mem_gb += 8 * nsph * model['max_layers'] # RCF
        host_mem_gb += 8 * (16 + 16 + 16 + 16) # CMULLR + CMUL + CEXTLR + CEXT
        host_mem_gb += 8 * 3 * 11 # UNIT VECTORS
        host_mem_gb += 8 + 8 # ARGI + ARGS
        host_mem_gb += 8 * (lm * 12 + nsph) # ZPV + GAPS
        host_mem_gb += 8 * (3 + 6 + 6) # GAP + GAPV + GAPM
        host_mem_gb += 16 * (6 + 6) # GAPP + GAPPM
        host_mem_gb += 8 * (2 * nsph + 4 * nsph) # TQSE + TQSPE + TQSS + TQSPS
        host_mem_gb += 8 * (2 * 4 + 4 * 4) # TQCE + TQCPE + TQCS + TQCPS
        host_mem_gb += 8 * (3 + 3) # TQEV + TQSV
        host_mem_gb += 16 * nsph * (half_layers + 1) # DC0M
        host_mem_gb += 8 * model['nxi'] # XIV
        host_mem_gb += 16 * 5 # ENT + ENTN + ARG + S0 + S0M
        
    host_mem_gb /= (1024.0 * 1024.0 * 1024.0)
    gpu_mem_gb /= (1024.0 * 1024.0 * 1024.0)
    print("INFO: 1 wavelength iteration uses %.5g Gb of host memory."%host_mem_gb)
    print("INFO: 1 wavelength iteration uses %.5g Gb of GPU memory"%gpu_mem_gb)
    print("      (or %.5g Gb, if iterative refinement is on)."%(3.0 * gpu_mem_gb))
    
## \brief Scan the calculation model.
#
#  The computational costs depennd the characteristics of the model. This
#  function scans the model configuration files and returns a model description
#  dictionary.
#
#  \returns model: `dict` A dictionary containing the model description.
def scan_model(edfb_name, geom_name):
    file_line = "INIT"
    model = {
        'app': "",
        'nsph': 0,
        'ies': 0,
        'li': 0,
        'le': 0,
        'npnt': 0,
        'npntts': 0,
        'nxi': 0,
        'configurations': 0,
        'max_layers': 0,
        'ros': [],
        'vec_nshl': [],
        'vec_types': [],
        'vec_x': [],
        'vec_y': [],
        'vec_z': []
    }
    # PARSING OF SCATTERER CONFIGURATION
    read_lines = 0
    max_layers = 0
    edfb_file = open(edfb_name, "r")
    file_line = edfb_file.readline()
    read_lines += 1
    iter_values = int_reg.finditer(file_line)
    array_values = iter_to_array(iter_values)
    if (len(array_values) != 2):
        print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
        print("   INVALID LINE: \"%s\""%file_line[:-1])
        print("   at line %d in %s."%(read_lines, edfb_name))
        edfb_file.close()
        return None
    model['nsph'] = int(array_values[0].group())
    model['ies'] = int(array_values[1].group())
    if (model['ies'] > 0):
        model['app'] = "INCLUSION"
        file_line = edfb_file.readline()
        read_lines += 1
        str_line = file_line.replace('D', 'E').replace('d', 'E')
        iter_values = number_reg.finditer(str_line)
        array_values = iter_to_array(iter_values)
        if (len(array_values) != 3):
            print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
            print("   INVALID LINE: \"%s\""%file_line[:-1])
            print("   at line %d in %s."%(read_lines, edfb_name))
            edfb_file.close()
            return None
        str_line = str_line[array_values[2].end():]
        iter_values = int_reg.finditer(str_line)
        array_values = iter_to_array(iter_values)
        if (len(array_values) != 4):
            print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
            print("   INVALID LINE: \"%s\""%file_line[:-1])
            print("   at line %d in %s."%(read_lines, edfb_name))
            edfb_file.close()
            return None
        idfc = int(array_values[0].group())
        nxi = int(array_values[1].group())
        instpc = int(array_values[2].group())
        insn = int(array_values[3].group())
        model['nxi'] = nxi
        if (instpc != 0):
            file_line.readline()
            file_line.readline()
            read_lines += 2
        else:
            for fi in range(nxi):
                edfb_file.readline()
            read_lines += nxi
        # end of if(instpc) block
        found_spheres = 0
        configurations = 0
        last_configuration = 0
        last_type = 0
        while (found_spheres < model['nsph']):
            file_line = edfb_file.readline()
            read_lines += 1
            iter_values = int_reg.finditer(file_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) < 1):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            found_spheres += len(array_values)
            for ci in range(len(array_values)):
                type_id = int(array_values[ci].group())
                if (type_id > last_type):
                    last_type = type_id
                    configurations += 1
                    last_configuration += 1
                model['vec_types'].append(last_configuration)
            # end for ci block
        # end while(found_spheres) block
        model['configurations'] = configurations
        for ci in range(configurations):
            file_line = edfb_file.readline()
            read_lines += 1
            iter_values = int_reg.finditer(file_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) < 1):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            nsh = int(array_values[0].group())
            model['vec_nshl'].append(nsh)
            if (ci == 0): nsh += 1
            if (nsh > max_layers): max_layers = nsh
            str_line = file_line[array_values[0].end():].replace('D', 'E').replace('d', 'E')
            iter_values = number_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 1):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            radius = float(array_values[0].group())
            for ish in range(nsh):
                file_line = edfb_file.readline()
                read_lines += 1
            if (ci == 0):
                str_line = file_line.replace('D', 'E').replace('d', 'E')
                iter_values = number_reg.finditer(str_line)
                array_values = iter_to_array(iter_values)
                if (len(array_values) != 1):
                    print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                    print("   INVALID LINE: \"%s\""%file_line[:-1])
                    print("   at line %d in %s."%(read_lines, edfb_name))
                    edfb_file.close()
                    return None
                factor = float(array_values[0].group())
                radius *= factor
            # end of if (ci == 0) block
            model['ros'].append(radius)
        # end of for ci block
        model['max_layers'] = max_layers
    else:
        # ies == 0
        if model['nsph'] == 1:
            model['app'] = "SPHERE"
            file_line = edfb_file.readline()
            str_line = file_line.replace('D', 'E').replace('d', 'E')
            iter_values = number_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 3):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            str_line = str_line[array_values[2].end():]
            iter_values = int_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 4):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            idfc = int(array_values[0].group())
            nxi = int(array_values[1].group())
            instpc = int(array_values[2].group())
            insn = int(array_values[3].group())
            model['nxi'] = nxi
            if (instpc != 0):
                file_line.readline()
                file_line.readline()
                read_lines += 2
            else:
                for fi in range(nxi):
                    edfb_file.readline()
                read_lines += nxi
            # end of if(instpc) block
            file_line = edfb_file.readline()
            file_line = edfb_file.readline()
            read_lines += 2
            model['vec_types'].append(0)
            iter_values = int_reg.finditer(file_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) < 1):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            max_layers = int(array_values[0].group())
            model['vec_nshl'].append(max_layers)
            str_line = file_line[array_values[0].end():].replace('D', 'E').replace('d', 'E')
            iter_values = number_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 1):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            radius = float(array_values[0].group())
            model['ros'].append(radius)
            model['max_layers'] = max_layers
        else:
            model['app'] = "CLUSTER"
            file_line = edfb_file.readline()
            read_lines += 1
            str_line = file_line.replace('D', 'E').replace('d', 'E')
            iter_values = number_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 3):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            str_line = str_line[array_values[2].end():]
            iter_values = int_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 4):
                print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                print("   INVALID LINE: \"%s\""%file_line[:-1])
                print("   at line %d in %s."%(read_lines, edfb_name))
                edfb_file.close()
                return None
            idfc = int(array_values[0].group())
            nxi = int(array_values[1].group())
            instpc = int(array_values[2].group())
            insn = int(array_values[3].group())
            model['nxi'] = nxi
            if (instpc != 0):
                file_line.readline()
                file_line.readline()
                read_lines += 2
            else:
                for fi in range(nxi):
                    edfb_file.readline()
                read_lines += nxi
            # end of if(instpc) block
            found_spheres = 0
            configurations = 0
            last_configuration = 0
            last_type = 0
            while (found_spheres < model['nsph']):
                file_line = edfb_file.readline()
                read_lines += 1
                iter_values = int_reg.finditer(file_line)
                array_values = iter_to_array(iter_values)
                if (len(array_values) < 1):
                    print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                    print("   INVALID LINE: \"%s\""%file_line[:-1])
                    print("   at line %d in %s."%(read_lines, edfb_name))
                    edfb_file.close()
                    return None
                found_spheres += len(array_values)
                for ci in range(len(array_values)):
                    type_id = int(array_values[ci].group())
                    if (type_id > last_type):
                        last_type = type_id
                        configurations += 1
                        last_configuration += 1
                    model['vec_types'].append(last_configuration)
                # end for ci block
            # end while(found_spheres) block
            model['configurations'] = configurations
            for ci in range(configurations):
                file_line = edfb_file.readline()
                read_lines += 1
                iter_values = int_reg.finditer(file_line)
                array_values = iter_to_array(iter_values)
                if (len(array_values) < 1):
                    print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                    print("   INVALID LINE: \"%s\""%file_line[:-1])
                    print("   at line %d in %s."%(read_lines, edfb_name))
                    edfb_file.close()
                    return None
                nsh = int(array_values[0].group())
                model['vec_nshl'].append(nsh)
                if (nsh > max_layers): max_layers = nsh
                str_line = file_line[array_values[0].end():].replace('D', 'E').replace('d', 'E')
                iter_values = number_reg.finditer(str_line)
                array_values = iter_to_array(iter_values)
                if (len(array_values) != 1):
                    print("ERROR: %s is not a valid scatterer configuration file!"%edfb_name)
                    print("   INVALID LINE: \"%s\""%file_line[:-1])
                    print("   at line %d in %s."%(read_lines, edfb_name))
                    edfb_file.close()
                    return None
                radius = float(array_values[0].group())
                for ish in range(nsh):
                    file_line = edfb_file.readline()
                    read_lines += 1
                model['ros'].append(radius)
            # end of for ci block
            model['max_layers'] = max_layers
        # end if model['nsph'] block
    # end if model['ies'] block
    edfb_file.close()
    
    # PARSING OF GEOMETRY CONFIGURATION
    read_lines = 0
    geom_file = open(geom_name, "r")
    file_line = geom_file.readline()
    read_lines += 1
    if (model['app'] != "SPHERE"):
        iter_values = int_reg.finditer(file_line)
        array_values = iter_to_array(iter_values)
        if (len(array_values) != 9):
            print("ERROR: %s is not a valid geometry configuration file!"%geom_name)
            print("   INVALID LINE: \"%s\""%file_line[:-1])
            print("   at line %d in %s."%(read_lines, geom_name))
            geom_file.close()
            return None
        nsph = int(array_values[0].group())
        if (nsph != model['nsph']):
            print("ERROR: %s is not consistent with %s!"%(geom_name, edfb_name))
            geom_file.close()
            return None
        model['li'] = int(array_values[1].group())
        model['le'] = int(array_values[2].group())
        model['npnt'] = int(array_values[5].group())
        model['npntts'] = int(array_values[6].group())
        for si in range(nsph):
            file_line = geom_file.readline()
            read_lines += 1
            str_line = file_line.replace('D', 'E').replace('d', 'E')
            iter_values = number_reg.finditer(str_line)
            array_values = iter_to_array(iter_values)
            if (len(array_values) != 3):
                print("ERROR: %s is not consistent with %s!"%(geom_name, edfb_name))
                geom_file.close()
                return None
            x = float(array_values[0].group())
            y = float(array_values[1].group())
            z = float(array_values[2].group())
            model['vec_x'].append(x)
            model['vec_y'].append(y)
            model['vec_z'].append(z)
    else:
        # model['app'] == "SPHERE"
        iter_values = int_reg.finditer(file_line)
        array_values = iter_to_array(iter_values)
        if (len(array_values) != 6):
            print("ERROR: %s is not a valid geometry configuration file!"%geom_name)
            print("   INVALID LINE: \"%s\""%file_line[:-1])
            print("   at line %d in %s."%(read_lines, geom_name))
            geom_file.close()
            return None
        nsph = int(array_values[0].group())
        if (nsph != model['nsph']):
            print("ERROR: %s is not consistent with %s!"%(geom_name, edfb_name))
            geom_file.close()
            return None
        model['li'] = int(array_values[1].group())
        model['npnt'] = int(array_values[3].group())
        model['npntts'] = int(array_values[4].group())
        model['vec_x'].append(0.0)
        model['vec_y'].append(0.0)
        model['vec_z'].append(0.0)
    geom_file.close()
    return model

## \brief Exit code (0 for success).
exit_code = main()
exit(exit_code)
