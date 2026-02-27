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

## @package model_maker
#  \brief Script to build models from YAML configuration files.
#
#  Script to assist in the creation of model input files starting
#  from a YAML descriptor.
#
#  The script requires python3.

import math
import multiprocessing
import numpy as np
import os
import pdb
import random
import yaml

## \brief 3D software generation capability flag.
allow_3d = True
try:
    import pyvista as pv
except ModuleNotFoundError as ex:
    print("WARNING: pyvista not found!")
    allow_3d = False

from pathlib import Path
from sys import argv

## \brief Main execution code
#
# `main()` is the function that handles the creation of the code configuration.
# It returns an integer value as exit code, using 0 to signal successful execution.
#
# \returns result: `int` Exit code (0 = SUCCESS).
def main():
    result = 0
    config = parse_arguments()
    if (config['help_mode'] or config['yml_file_name'] == ""):
        print_help()
    else:
        rnd_attempts = config['attempt_limit']
        sconf, gconf = load_model(config['yml_file_name'], rnd_attempts)
        if (sconf is not None) and (gconf is not None):
            result = write_legacy_sconf(sconf)
            if (result == 0):
                result += write_legacy_gconf(gconf)
            if (result == 0):
                print_model_summary(sconf, gconf)
        else:
            print("ERROR: could not create configuration.")
            result = 1
    return result

## \brief Populate the dielectric constant data via interpolation
#
#  \param sconf: `dict` Scatterer configuration dictionary.
#  \return result: `int` An exit code (0 if successful).
def interpolate_constants(sconf):
    result = 0
    err_arg = ""
    try:
        for i in range(sconf['configurations']):
            for j in range(sconf['nshl'][i]):
                file_idx = sconf['dielec_id'][i][j]
                dielec_path = Path(sconf['dielec_path'], sconf['dielec_file'][int(file_idx) - 1])
                err_arg = str(dielec_path)
                file_name = err_arg
                dielec_file = open(file_name, 'r')
                wavelengths = []
                rpart = []
                ipart = []
                str_line = dielec_file.readline()
                while (str_line != ""):
                    if (not str_line.startswith('#')):
                        split_line = str_line.split(',')
                        if (len(split_line) == 3):
                            wavelengths.append(float(split_line[0]))
                            rpart.append(float(split_line[1]))
                            ipart.append(float(split_line[2]))
                    str_line = dielec_file.readline()
                dielec_file.close()
                wi = 0
                x0 = 0.0
                x1 = 0.0
                ry0 = 0.0
                iy0 = 0.0
                ry1 = 0.0
                iy1 = 0.0
                for dci in range(sconf['nxi']):
                    w = sconf['vec_xi'][dci]
                    while (w > x1):
                        x0 = wavelengths[wi]
                        ry0 = rpart[wi]
                        iy0 = ipart[wi]
                        if (wi == len(wavelengths)):
                            print("ERROR: file %s does not cover requested wavelengths!"%file_name)
                            return 1
                        wi += 1
                        x1 = wavelengths[wi]
                        ry1 = rpart[wi]
                        iy1 = ipart[wi]
                    if (wi > 0):
                        x0 = wavelengths[wi - 1]
                        ry0 = rpart[wi - 1]
                        iy0 = ipart[wi - 1]
                        dx = w - x0
                        dry = (ry1 - ry0) / (x1 - x0) * dx
                        diy = (iy1 - iy0) / (x1 - x0) * dx
                        ry = ry0 + dry
                        iy = iy0 + diy
                        sconf['rdc0'][j][i][dci] = ry
                        sconf['idc0'][j][i][dci] = iy
                    else:
                        if (wavelengths[wi] == w):
                            sconf['rdc0'][j][i][dci] = rpart[0]
                            sconf['idc0'][j][i][dci] = ipart[0]
                        else:
                            print("ERROR: file %s does not cover requested wavelengths!"%file_name)
                            return 2
    except FileNotFoundError as ex:
        print("ERROR: file not found %s!"%err_arg)
        return 3
    return result

## \brief Create tha calculation configuration structure from YAML input.
#
#  \param model_file: `str` Full path to the YAML input file.
#  \param rnd_attempts: `int` Limiting number of attempts to randomly place a sphere.
#  \return sconf, gconf: `tuple` A dictionary tuple for scatterer and
#                        geometric configurations.
def load_model(model_file, rnd_attempts):
    sconf = None
    gconf = None
    model = None
    try:
        with open(model_file, 'r') as stream:
            model = yaml.safe_load(stream)
    except yaml.YAMLError:
        print("ERROR: " + model_file + " is not a valid YAML file!")
    except FileNotFoundError:
        print("ERROR: " + model_file + " was not found!")
    if model is not None:
        max_rad = 0.0
        make_3d = False
        try:
            make_3d = False if model['system_settings']['make_3D'] == "0" else True
        except KeyError:
            make_3d = False
        if (make_3d and not allow_3d):
            print("WARNING: 3D visualization of models is not available. Disabling.")
            make_3d = False
        # Create the sconf dict
        sconf = {
            'out_file': Path(
                model['input_settings']['input_folder'],
                model['input_settings']['spheres_file']
            )
        }
        sconf['nsph'] = int(model['particle_settings']['n_spheres'])
        sconf['application'] = model['particle_settings']['application']
        sconf['ies'] = 1 if sconf['application'] == "INCLUSION" else 0
        sconf['exdc'] = float(model['material_settings']['extern_diel'])
        sconf['wp'] = float(model['radiation_settings']['wp'])
        sconf['xip'] = float(model['radiation_settings']['xip'])
        sconf['idfc'] = int(model['material_settings']['diel_flag'])
        sconf['instpc'] = int(model['radiation_settings']['step_flag'])
        sconf['xi_start'] = float(model['radiation_settings']['scale_start'])
        sconf['xi_end'] = float(model['radiation_settings']['scale_end'])
        sconf['xi_step'] = float(model['radiation_settings']['scale_step'])
        sconf['configurations'] = int(model['particle_settings']['n_types'])
        sconf['dielec_path'] = model['material_settings']['dielec_path']
        sconf['dielec_file'] = model['material_settings']['dielec_file']
        num_dielec = 0 # len(model['particle_settings']['dielec_id'])
        if (sconf['idfc'] == -1):
            num_dielec = len(model['material_settings']['diel_func'])
        elif (sconf['idfc'] == 0):
            num_dielec = len(model['particle_settings']['dielec_id'])
        if (len(model['particle_settings']['n_layers']) != sconf['configurations']):
            print("ERROR: Declared number of layers does not match number of types!")
            return (None, None)
        else:
            sconf['nshl'] = [0 for i in range(sconf['configurations'])]
            for i in range(sconf['configurations']):
                sconf['nshl'][i] = int(model['particle_settings']['n_layers'][i])
        max_layers = max(sconf['nshl'])
        if (sconf['application'] == "INCLUSION"):
            if (max_layers < sconf['nshl'][0] + 1):
                max_layers = sconf['nshl'][0] + 1
        if (num_dielec != sconf['configurations']):
            print("ERROR: declared array of optical constants does not match configurations!")
            return (None, None)
        else:
            if (sconf['idfc'] == 0):
                sconf['dielec_id'] = [
                    [ 0 for j in range(max_layers)] for i in range(sconf['configurations'])
                ]
                for i in range(sconf['configurations']):
                    expected_layer_num = 1 + int(sconf['nshl'][i] / 2)
                    if (sconf['application'] == "INCLUSION" and i == 0):
                        expected_layer_num += 1
                    if (len(model['particle_settings']['dielec_id'][i]) != expected_layer_num):
                        print("ERROR: Declared materials in type %d do not match the number of layers!"%(i + 1))
                        return (None, None)
                    else:
                        for j in range(1 + int(sconf['nshl'][i] / 2)):
                            sconf['dielec_id'][i][j] = float(model['particle_settings']['dielec_id'][i][j])
        if (model['radiation_settings']['scale_name'] == "XI"):
            sconf['insn'] = 1
            sconf['nxi'] = 1 + int((sconf['xi_end'] - sconf['xi_start']) / sconf['xi_step'])
            sconf['vec_xi'] = [0.0 for i in range(sconf['nxi'])]
            for i in range(sconf['nxi']):
                sconf['vec_xi'][i] = sconf['xi_start'] + i * sconf['xi_step']
            # Set up dielectric constants
            allow_dielec = True # TODO: define logic to check dielectric constants
            if (not allow_dielec):
                print("ERROR: delcared dielectric constants do not match number of sphere types!")
                return (None, None)
            else:
                if (sconf['idfc'] == -1):
                    sconf['rdc0'] = [
                        [
                            [0.0 for k in range(1)] for j in range(sconf['configurations'])
                        ] for i in range(max_layers)
                    ]
                    sconf['idc0'] = [
                        [
                            [0.0 for k in range(1)] for j in range(sconf['configurations'])
                        ] for i in range(max_layers)
                    ]
                    for di in range(max_layers):
                        for dj in range(sconf['configurations']):
                            if (len(model['material_settings']['diel_func'][dj]) / 2 != sconf['nshl'][dj]):
                                print("ERROR: dielectric constants for type %d do not match number of layers!"%(dj + 1))
                                return (None, None)
                            else:
                                sconf['rdc0'][di][dj][0] = float(model['material_settings']['diel_func'][dj][2 * di])
                                sconf['idc0'][di][dj][0] = float(model['material_settings']['diel_func'][dj][2 * di + 1])
                else: # sconf[idfc] != -1 and scaling on XI
                    print("ERROR: for scaling on XI, optical constants must be defined!")
                    return (None, None)
        elif (model['radiation_settings']['scale_name'] == "WAVELENGTH"):
            sconf['insn'] = 3
            if (model['material_settings']['match_mode'] == "INTERPOLATE"):
                sconf['nxi'] = 1 + int((sconf['xi_end'] - sconf['xi_start']) / sconf['xi_step'])
                sconf['vec_xi'] = [0.0 for i in range(sconf['nxi'])]
                for i in range(sconf['nxi']):
                    sconf['vec_xi'][i] = sconf['xi_start'] + i * sconf['xi_step']
                # Set up the dielectric constants
                if (sconf['idfc'] == 0):
                    sconf['rdc0'] = [
                        [
                            [0.0 for k in range(sconf['nxi'])] for j in range(sconf['configurations'])
                        ] for i in range(max_layers)
                    ]
                    sconf['idc0'] = [
                        [
                            [0.0 for k in range(sconf['nxi'])] for j in range(sconf['configurations'])
                        ] for i in range(max_layers)
                    ]
                    check = interpolate_constants(sconf)
                    if (check != 0):
                        return (None, None)
                else: # sconf[idfc] != 0 and scaling on wavelength
                    print("ERROR: for wavelength scaling, optical constants must be tabulated!")
                    return (None, None)
            elif (model['material_settings']['match_mode'] == "GRID"):
                check = match_grid(sconf)
                if (check != 0):
                    return(None, None)
            else:
                print("ERROR: %s is not a recognized match mode!"%(model['material_settings']['match_mode']))
                return (None, None)
        else:
            print("ERROR: %s is not a supported scaling mode!"%(model['radiation_settings']['scale_name']))
            return (None, None)
        sph_types = model['particle_settings']['sph_types']
        if (len(sph_types) == sconf['nsph']):
            sconf['vec_types'] = [int(str_typ) for str_typ in sph_types]
        else:
            if (len(sph_types) != 0):
                print("ERROR: vector of sphere types does not match the declared number of spheres!")
                return (None, None)
            else: # len(sph_types) = 0
                len_vec_x = len(model['geometry_settings']['x_coords'])
                len_vec_y = len(model['geometry_settings']['y_coords'])
                len_vec_z = len(model['geometry_settings']['z_coords'])
                if (len_vec_x != 0 or len_vec_y != 0 or len_vec_z != 0):
                    print("ERROR: cannot assign random types with explicit sphere positions!")
                    return (None, None)
                else:
                    sconf['vec_types'] = [0 for ti in range(sconf['nsph'])]
        if (len(model['particle_settings']['radii']) != sconf['configurations']):
            print("ERROR: Declared number of radii does not match number of types!")
            return (None, None)
        else:
            sconf['ros'] = [0.0 for i in range(sconf['configurations'])]
            for i in range(sconf['configurations']):
                sconf['ros'][i] = float(model['particle_settings']['radii'][i])
        if (len(model['particle_settings']['rad_frac']) != sconf['configurations']):
            print("ERROR: Declared number of fractional radii does not match number of types!")
            return (None, None)
        else:
            sconf['rcf'] = [
                [0.0 for j in range(max_layers)] for i in range(sconf['configurations'])
            ]
            for i in range(sconf['configurations']):
                expected_radii = sconf['nshl'][i]
                if (sconf['application'] == "INCLUSION" and i == 0):
                    expected_radii += 1
                if (len(model['particle_settings']['rad_frac'][i]) != expected_radii):
                    print("ERROR: Declared transition radii in type %d do not match the number of layers!"%(i + 1))
                    return (None, None)
                else:
                    expected_radii = sconf['nshl'][i]
                    if (sconf['application'] == "INCLUSION" and i == 0):
                        expected_radii += 1
                    for j in range(expected_radii):
                        sconf['rcf'][i][j] = float(model['particle_settings']['rad_frac'][i][j])
        # Create the gconf dict
        use_refinement = True
        dyn_orders = True
        try:
            use_refinement = False if int(model['runtime']['refinement']) == 0 else True
        except KeyError:
            use_refinement = True
        try:
            dyn_orders = False if int(model['runtime']['dyn_orders']) == 0 else True
        except KeyError:
            dyn_orders = True
        str_polar = model['radiation_settings']['polarization']
        if (str_polar not in ["LINEAR", "CIRCULAR"]):
            print("ERROR: %s is not a recognized polarization state."%str_polar)
            return (None, None)
        gconf = {
            'out_file': Path(
                model['input_settings']['input_folder'],
                model['input_settings']['geometry_file']
            )
        }
        gconf['use_refinement'] = use_refinement
        gconf['dyn_orders'] = dyn_orders
        gconf['max_host_ram_gb'] = 0
        try:
            gconf['max_host_ram_gb'] = float(model['system_settings']['max_host_ram'])
        except KeyError as ex:
            print("WARNING: no host RAM declared. Cannot estimate recommended execution.")
        gconf['max_gpu_ram_gb'] = 0
        try:
            gconf['max_gpu_ram_gb'] = float(model['system_settings']['max_gpu_ram'])
        except KeyError as ex:
            print("WARNING: no GPU RAM declared.")
        gconf['nsph'] = sconf['nsph']
        gconf['application'] = model['particle_settings']['application']
        gconf['li'] = int(model['geometry_settings']['li'])
        gconf['le'] = int(
            gconf['li'] if gconf['application'] == "SPHERE" else model['geometry_settings']['le']
        )
        gconf['inpol'] = 0 if str_polar == "LINEAR" else 1
        gconf['npnt'] = int(model['geometry_settings']['npnt'])
        gconf['npntts'] = int(model['geometry_settings']['npntts'])
        if (gconf['application'] != "SPHERE"):
            gconf['iavm'] = int(model['geometry_settings']['iavm'])
        gconf['isam'] = int(model['geometry_settings']['isam'])
        gconf['jwtm'] = int(model['output_settings']['jwtm'])
        gconf['th'] = float(model['geometry_settings']['in_th_start'])
        gconf['thstp'] = float(model['geometry_settings']['in_th_step'])
        gconf['thlst'] = float(model['geometry_settings']['in_th_end'])
        gconf['ph'] = float(model['geometry_settings']['in_ph_start'])
        gconf['phstp'] = float(model['geometry_settings']['in_ph_step'])
        gconf['phlst'] = float(model['geometry_settings']['in_ph_end'])
        gconf['ths'] = float(model['geometry_settings']['sc_th_start'])
        gconf['thsstp'] = float(model['geometry_settings']['sc_th_step'])
        gconf['thslst'] = float(model['geometry_settings']['sc_th_end'])
        gconf['phs'] = float(model['geometry_settings']['sc_ph_start'])
        gconf['phsstp'] = float(model['geometry_settings']['sc_ph_step'])
        gconf['phslst'] = float(model['geometry_settings']['sc_ph_end'])
        gconf['vec_sph_x'] = [0.0 for i in range(gconf['nsph'])]
        gconf['vec_sph_y'] = [0.0 for i in range(gconf['nsph'])]
        gconf['vec_sph_z'] = [0.0 for i in range(gconf['nsph'])]
        if (gconf['application'] != "SPHERE" or gconf['nsph'] != 1):
            len_vec_x = len(model['geometry_settings']['x_coords'])
            len_vec_y = len(model['geometry_settings']['y_coords'])
            len_vec_z = len(model['geometry_settings']['z_coords'])
            if (len_vec_x != len_vec_y):
                print("ERROR: X and Y coordinate vectors have different lengths!")
                return (None, None)
            if (len_vec_x != len_vec_z):
                print("ERROR: X and Z coordinate vectors have different lengths!")
                return (None, None)
            if (len_vec_y != len_vec_z):
                print("ERROR: Y and Z coordinate vectors have different lengths!")
                return (None, None)
            if (len_vec_x == 0):
                # Generate random cluster
                rnd_seed = int(model['system_settings']['rnd_seed'])
                try:
                    max_rad = float(model['particle_settings']['max_rad'])
                except KeyError:
                    print("ERROR: random model generation requires specification of particle_settings:max_rad.")
                    return (None, None)
                rnd_engine = "COMPACT"
                try:
                    rnd_engine = model['system_settings']['rnd_engine']
                except KeyError:
                    # use compact generator, if no specification is given
                    rnd_engine = "COMPACT"
                if (rnd_engine == "COMPACT"):
                    check = random_compact(sconf, gconf, rnd_seed, max_rad, rnd_attempts)
                    if (check == -1):
                        print("ERROR: compact random generator works only when all sphere types have the same radius.")
                        return (None, None)
                    elif (check == -2):
                        print("ERROR: sub-particle radius larger than particle radius.")
                        return (None, None)
                    elif (check == -3):
                        print("ERROR: requested number of spheres cannot fit in allowed volume.")
                        return (None, None)
                elif (rnd_engine == "LOOSE"):
                    # random_aggregate() checks internally whether application is INCLUSION
                    check = random_aggregate(sconf, gconf, rnd_seed, max_rad, rnd_attempts)
                else:
                    print("ERROR: unrecognized random generator engine.")
                    return (None, None)
                if (check != sconf['nsph']):
                    print("WARNING: placed only %d out of %d requested spheres."%(check, sconf['nsph']))
                    sconf['nsph'] = check
                    gconf['nsph'] = check
            else:
                if (len(model['geometry_settings']['x_coords']) != gconf['nsph']):
                    print("ERROR: coordinate vectors do not match the number of spheres!")
                    return (None, None)
                for si in range(gconf['nsph']):
                    gconf['vec_sph_x'][si] = float(model['geometry_settings']['x_coords'][si])
                    gconf['vec_sph_y'][si] = float(model['geometry_settings']['y_coords'][si])
                    gconf['vec_sph_z'][si] = float(model['geometry_settings']['z_coords'][si])
        if (make_3d and allow_3d):
            if (max_rad == 0.0):
                max_rad = 20.0 * max(sconf['ros'])
            write_obj(sconf, gconf, max_rad)
        test_system_resources(model, gconf, sconf)
    else: # model is None
        print("ERROR: could not parse " + model_file + "!")
    return (sconf, gconf)

## \brief Populate the dielectric constant data matching a grid.
#
#  Important note: if the configuration requests that more than one
#  optical constants file should be used, all the files must provide
#  their constants for the same vector of wavelengths.
#
#  \param sconf: `dict` Scatterer configuration dictionary.
#  \return result: `int` An exit code (0 if successful).
def match_grid(sconf):
    result = 0
    max_layers = 0
    nxi = 0
    sconf['vec_xi'] = []
    err_arg = ""
    try:
        for i in range(sconf['configurations']):
            layers = sconf['nshl'][i]
            if (sconf['application'] == "INCLUSION" and i == 0):
                layers += 1
            for j in range(layers):
                file_idx = sconf['dielec_id'][i][j]
                dielec_path = Path(sconf['dielec_path'], sconf['dielec_file'][int(file_idx) - 1])
                err_arg = str(dielec_path)
                file_name = err_arg
                dielec_file = open(file_name, 'r')
                wavelengths = []
                rpart = []
                ipart = []
                str_line = dielec_file.readline()
                while (str_line != ""):
                    if (not str_line.startswith('#')):
                        split_line = str_line.split(',')
                        if (len(split_line) == 3):
                            wavelengths.append(float(split_line[0]))
                            rpart.append(float(split_line[1]))
                            ipart.append(float(split_line[2]))
                    str_line = dielec_file.readline()
                dielec_file.close()
                if (max_layers == 0):
                    # This is executed only once
                    max_layers = max(sconf['nshl'])
                    if (sconf['application'] == "INCLUSION" and max_layers < sconf['nshl'][0] + 1):
                        max_layers = sconf['nshl'][0] + 1
                    w_start = sconf['xi_start']
                    w_end = sconf['xi_end']
                    for wi in range(len(wavelengths)):
                        w = wavelengths[wi]
                        if (w >= w_start and w <= w_end):
                            sconf['vec_xi'].append(w)
                            nxi += 1
                    sconf['rdc0'] = [
                        [
                            [
                                0.0 for dk in range(nxi)
                            ] for dj in range(sconf['configurations'])
                        ] for di in range(max_layers)
                    ]
                    sconf['idc0'] = [
                        [
                            [
                                0.0 for dk in range(nxi)
                            ] for dj in range(sconf['configurations'])
                        ] for di in range(max_layers)
                    ]
                    sconf['nxi'] = nxi
                # This is executed for all layers in all configurations
                wi = 0
                x = wavelengths[wi]
                ry = rpart[wi]
                iy = ipart[wi]
                for dci in range(sconf['nxi']):
                    w = sconf['vec_xi'][dci]
                    while (w > x):
                        x = wavelengths[wi]
                        ry = rpart[wi]
                        iy = ipart[wi]
                        if (wi == len(wavelengths)):
                            print("ERROR: file %s does not cover requested wavelengths!"%file_name)
                            return 1
                        wi += 1
                    sconf['rdc0'][j][i][dci] = ry
                    sconf['idc0'][j][i][dci] = iy
    except FileNotFoundError as ex:
        print("ERROR: file not found %s!"%err_arg)
        return 3
    return result

## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of optional arguments.
#  The purpose of this function is to parse the command line in search for
#  such arguments and prepare the execution accordingly.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'attempt_limit': 100,
        'yml_file_name': "",
        'help_mode': False,
    }
    yml_set = False
    for arg in argv[1:]:
        if (arg.startswith("--help")):
            config['help_mode'] = True
        elif (arg.startswith("--attempts=")):
            split_arg = arg.split('=')
            config['attempt_limit'] = int(split_arg[1])
        elif (not yml_set):
            if (not arg.startswith("--")):
                config['yml_file_name'] = arg
                yml_set = True
        else:
            raise Exception("Unrecognized argument \'{0:s}\'".format(arg))
    return config

## \brief Print a command-line help summary.
def print_help():
    print("###############################################           ")
    print("#                                             #           ")
    print("#           NPtm_code MODEL_MAKER             #           ")
    print("#                                             #           ")
    print("###############################################           ")
    print("                                                          ")
    print("Create input files for FORTRAN and C++ code.              ")
    print("                                                          ")
    print("Usage: \"./model_maker.py CONFIG [OPTIONS]\"              ")
    print("                                                          ")
    print("CONFIG must be a valid YAML configuration file.           ")
    print("                                                          ")
    print("Valid options are:                                        ")
    print("--attempts=N          Try to place a sphere up to N times.")
    print("--help                Print this help and exit.           ")
    print("                                                          ")

## \brief Print a summary of model properties.
#
#  This function provides a summary of useful information concerning the
#  radii of the particle monomers and of the equivalent mass sphere, to
#  assist in the selection of the proper starting orders.
#
#  \param scatterer: `dict` A dictionary for the scatterer configuration.
#  \param geometry: `dict` A dictionary for the geometry configuration.
def print_model_summary(scatterer, geometry):
    avgX = 0.0
    avgY = 0.0
    avgZ = 0.0
    Rmin = 0.0
    Rmax = 0.0
    Reqm = 0.0
    R3tot = 0.0
    Rcirc = 0.0
    square_farthest = 0.0
    for i in range(scatterer['nsph']):
        avgX += geometry['vec_sph_x'][i]
        avgY += geometry['vec_sph_y'][i]
        avgZ += geometry['vec_sph_z'][i]
        sph_type_index = scatterer['vec_types'][i] - 1
        ros = scatterer['ros'][sph_type_index]
        R3tot += math.pow(ros, 3.0)
        if (ros > Rmax):
            Rmax = ros
        if (Rmin == 0.0 or ros < Rmin):
            Rmin = ros
    Reqm = math.pow(R3tot, 1.0 / 3.0)
    avgX /= scatterer['nsph']
    avgY /= scatterer['nsph']
    avgZ /= scatterer['nsph']
    for i in range(scatterer['nsph']):
        sph_type_index = scatterer['vec_types'][i] - 1
        ros = scatterer['ros'][sph_type_index]
        dX = geometry['vec_sph_x'][i] - avgX
        dY = geometry['vec_sph_y'][i] - avgY
        dZ = geometry['vec_sph_z'][i] - avgZ
        square_range = dX * dX + dY * dY + dZ * dZ + ros * ros
        if (square_range > square_farthest):
            square_farthest = square_range
    Rcirc = math.sqrt(square_farthest)
    if (Rmax == Rmin):
        print("INFO: monomer radius is Rsph = %.5em"%Rmin)
    else:
        print("INFO: smallest monomer radius is Rmin = %.5em"%Rmin)
        print("INFO: largest monomer radius is Rmax = %.5em"%Rmax)
    print("INFO: equivalent volume radius is Reqv = %.5em"%Reqm)
    print("INFO: minimum encircling radius is Rcirc = %.5em"%Rcirc)

## \brief Generate a random cluster aggregate from YAML configuration options.
#
#  This function generates a random aggregate of spheres using radial ejection
#  in random directions of new spheres until they become tangent to the
#  outermost sphere existing in that direction. The result of the generated
#  model is directly saved in the parameters of the scatterer and geometry
#  configuration dictionaries.
#
#  \param scatterer: `dict` Scatterer configuration dictionary (gets modified)
#  \param geometry: `dict` Geometry configuration dictionary (gets modified)
#  \param seed: `int` Seed for the random sequence generation
#  \param max_rad: `float` Maximum allowed radial extension of the aggregate
#  \param max_attempts: `int` Maximum number of attempts to place a particle in any direction
#  \return result: `int` Function exit code (0 for success, otherwise number of
#  spheres that could not be placed)
def random_aggregate(scatterer, geometry, seed, max_rad, max_attempts=100):
    result = 0
    random.seed(seed)
    nsph = scatterer['nsph']
    vec_thetas = [0.0 for i in range(nsph)]
    vec_phis = [0.0 for i in range(nsph)]
    vec_rads = [0.0 for i in range(nsph)]
    vec_types = []
    n_types = scatterer['configurations']
    if (0 in scatterer['vec_types']):
        tincrement = 1 if scatterer['application'] != "INCLUSION" else 2
        for ti in range(nsph):
            itype = tincrement + int(n_types * random.random())
            scatterer['vec_types'][ti] = itype
        if (scatterer['application'] == "INCLUSION"):
            scatterer['vec_types'][0] = 1
    sph_type_index = scatterer['vec_types'][0] - 1
    vec_spheres = [{'itype': sph_type_index + 1, 'x': 0.0, 'y': 0.0, 'z': 0.0}]
    vec_rads[0] = scatterer['ros'][sph_type_index]
    vec_types.append(sph_type_index + 1)
    placed_spheres = 1
    attempts = 0
    for i in range(1, nsph):
        sph_type_index = scatterer['vec_types'][i] - 1
        vec_rads[i] = scatterer['ros'][sph_type_index]
        is_placed = False
        while (not is_placed):
            if (attempts > max_attempts):
                result += 1
                break # while(not is_placed)
            vec_thetas[i] = math.pi * random.random()
            vec_phis[i] = 2.0 * math.pi * random.random()
            rho = vec_rads[0] + vec_rads[i]
            z = rho * math.cos(vec_thetas[i])
            y = rho * math.sin(vec_thetas[i]) * math.sin(vec_phis[i])
            x = rho * math.sin(vec_thetas[i]) * math.cos(vec_phis[i])
            j = 0
            while (j < i - 1):
                j += 1
                dx2 = (x - vec_spheres[j]['x']) * (x - vec_spheres[j]['x'])
                dy2 = (y - vec_spheres[j]['y']) * (y - vec_spheres[j]['y'])
                dz2 = (z - vec_spheres[j]['z']) * (z - vec_spheres[j]['z'])
                dist2 = dx2 + dy2 + dz2
                rr2 = (vec_rads[i] + vec_rads[j]) * (vec_rads[i] + vec_rads[j])
                if (dist2 < 0.9999 * rr2):
                    # Spheres i and j are compenetrating.
                    # Sphere i is moved out radially until it becomes externally
                    # tangent to sphere j. Then the check is repeated, to verify
                    # that no other sphere was penetrated. The process is iterated
                    # until sphere i is placed or the maximum allowed radius is
                    # reached.
                    # breakpoint()
                    sinthi = math.sin(vec_thetas[i])
                    sinthj = math.sin(vec_thetas[j])
                    costhi = math.cos(vec_thetas[i])
                    costhj = math.cos(vec_thetas[j])
                    sinphi = math.sin(vec_phis[i])
                    sinphj = math.sin(vec_phis[j])
                    cosphi = math.cos(vec_phis[i])
                    cosphj = math.cos(vec_phis[j])
                    cosalpha = (
                        sinthi * cosphi * sinthj * cosphj
                        + sinthi * sinphi * sinthj * sinphj
                        + costhi * costhj
                    )
                    D12 = math.sqrt(
                        vec_spheres[j]['x'] * vec_spheres[j]['x']
                        + vec_spheres[j]['y'] * vec_spheres[j]['y']
                        + vec_spheres[j]['z'] * vec_spheres[j]['z']
                    )
                    O1K = D12 * cosalpha
                    sinalpha = math.sqrt(1.0 - cosalpha * cosalpha)
                    sinbetaprime = D12 / (vec_rads[i] + vec_rads[j]) * sinalpha
                    cosbetaprime = math.sqrt(1.0 - sinbetaprime * sinbetaprime)
                    Op3K = (vec_rads[i] + vec_rads[j]) * cosbetaprime
                    rho = O1K + Op3K
                    z = rho * math.cos(vec_thetas[i])
                    y = rho * math.sin(vec_thetas[i]) * math.sin(vec_phis[i])
                    x = rho * math.sin(vec_thetas[i]) * math.cos(vec_phis[i])
                    j = 0
                    continue # while(j < i - 1)
            if (rho + vec_rads[i] > max_rad):
                # The current direction is filled. Try another one.
                attempts += 1
                continue # while(not is_placed)
            vec_spheres.append({
                'itype': sph_type_index + 1,
                'x': x,
                'y': y,
                'z': z
            })
            is_placed = True
            placed_spheres += 1
            vec_types.append(sph_type_index + 1)
            attempts = 0
        # end while loop
    # end for i loop
    scatterer['vec_types'] = vec_types
    sph_index = 0
    for sphere in sorted(vec_spheres, key=lambda item: item['itype']):
        scatterer['vec_types'][sph_index] = sphere['itype']
        geometry['vec_sph_x'][sph_index] = sphere['x']
        geometry['vec_sph_y'][sph_index] = sphere['y']
        geometry['vec_sph_z'][sph_index] = sphere['z']
        sph_index += 1
    result = placed_spheres
    return result

## \brief Generate a random compact cluster from YAML configuration options.
#
#  This function generates a random aggregate of spheres using the maximum
#  compactness packaging to fill a spherical volume with given maximum radius,
#  then it proceeds by subtracting random spheres from the outer layers, until
#  the aggregate is reduced to the desired number of spheres. The function
#  can only be used if all sphere types have the same radius. The result of the
#  generated model is directly saved in the parameters of the scatterer and
#  geometry configuration dictionaries.
#
#  \param scatterer: `dict` Scatterer configuration dictionary (gets modified)
#  \param geometry: `dict` Geometry configuration dictionary (gets modified)
#  \param seed: `int` Seed for the random sequence generation
#  \param max_rad: `float` Maximum allowed radial extension of the aggregate
#  \return result: `int` Function exit code (0 for success, otherwise error code)
def random_compact(scatterer, geometry, seed, max_rad):
    result = 0
    random.seed(seed)
    nsph = scatterer['nsph']
    n_types = scatterer['configurations']
    radius = scatterer['ros'][0]
    # Return an error code if types have different radii
    if (max(scatterer['ros']) != min(scatterer['ros'])):
        result = -1
    elif (radius > max_rad):
        # Requested spheres are larger than the maximum allowed volume.
        # End function with error code -2.
        result = -2
    else:
        x_centers = np.arange(-1.0 * max_rad + 2.0 * radius, max_rad, 2.0 * radius)
        x_size = len(x_centers)
        y_size = int(2.0 * max_rad / ((1.0 + math.sqrt(3.0) / 3.0) * radius))
        z_size = int(2.0 * max_rad / ((1.0 + 2.0 * math.sqrt(6.0) / 3.0) * radius))
        tmp_spheres = []
        n_cells = x_size * y_size * z_size
        print("INFO: the cubic space would contain %d spheres."%n_cells)
        k = 0
        z = -max_rad + radius
        while (z < max_rad - radius):
            j = 0
            y = -max_rad + radius
            while (y < max_rad - radius):
                for i in range(len(x_centers)):
                    x = (2 * (i + 1) + (j + k) % 2) * radius - max_rad
                    extent = radius + math.sqrt(x * x + y * y + z * z)
                    if (extent < max_rad):
                        tmp_spheres.append({
                            'itype': 1,
                            'x': x,
                            'y': y,
                            'z': z
                        })
                #
                j += 1
                y = math.sqrt(3.0) * (j + (k % 2) / 3.0) * radius - max_rad + radius
            k += 1
            z = 2.0 / 3.0 * math.sqrt(6.0) * k * radius - max_rad + radius
        #tmp_spheres = [{'itype': 1, 'x': 0.0, 'y': 0.0, 'z': 0.0}]
        current_n = len(tmp_spheres)
        print("INFO: before erosion there are %d spheres in use."%current_n)
        rho = 10.0 * max_rad
        discard_rad = 100.0 * max_rad
        while (current_n > nsph):
            theta = math.pi * random.random()
            phi = 2.0 * math.pi * random.random()
            x0 = rho * math.sin(theta) * math.cos(phi)
            y0 = rho * math.sin(theta) * math.sin(phi)
            z0 = rho * math.cos(theta)
            closest_index = 0
            minimum_distance = 1000.0 * max_rad
            for di in range(len(tmp_spheres)):
                x1 = tmp_spheres[di]['x']
                if (x1 == discard_rad):
                    continue
                y1 = tmp_spheres[di]['y']
                z1 = tmp_spheres[di]['z']
                distance = math.sqrt(
                    (x1 - x0) * (x1 - x0)
                    + (y1 - y0) * (y1 - y0)
                    + (z1 - z0) * (z1 - z0)
                )
                if (distance < minimum_distance):
                    closest_index = di
                    minimum_distance = distance
            tmp_spheres[closest_index]['x'] = discard_rad
            current_n -= 1
        vec_spheres = []
        sph_index = 0
        # Generate a vector of types if none is given
        if (0 in scatterer['vec_types']):
            tincrement = 1 if scatterer['application'] != "INCLUSION" else 2
            for ti in range(current_n):
                itype = tincrement + int(n_types * random.random())
                scatterer['vec_types'][ti] = itype
            if (scatterer['application'] == "INCLUSION"):
                scatterer['vec_types'][0] = 1
        for ti in range(len(tmp_spheres)):
            sphere = tmp_spheres[ti]
            if (sphere['x'] < max_rad):
                sphere['itype'] = scatterer['vec_types'][sph_index]
                sph_index += 1
                vec_spheres.append(sphere)
        sph_index = 0
        for sphere in sorted(vec_spheres, key=lambda item: item['itype']):
            scatterer['vec_types'][sph_index] = sphere['itype']
            geometry['vec_sph_x'][sph_index] = sphere['x']
            geometry['vec_sph_y'][sph_index] = sphere['y']
            geometry['vec_sph_z'][sph_index] = sphere['z']
            sph_index += 1
    return current_n

## \brief Perform a preliminary test of the resources required by the model.
#
#  The resolution of models requires the availability of memory resources
#  that increase as a function of the model complexity. This function aims
#  at evaluating the complexity of the model and estimate whether the
#  declared system resources are sufficient to run the calculation.
#
#  \param model: `dict` Model description dictionary.
#  \param gconf: `dict` Geometry description dictionary.
#  \param sconf: `dict` Scattering description dictionary.
def test_system_resources(model, gconf, sconf):
    le = gconf['le'] if (gconf['application'] != "SPH") else 0
    li = gconf['li']
    lm = le if le > li else li
    nsph = gconf['nsph']
    nlim = li * (li + 2)
    nlem = le * (le + 2)
    nlemt = 2 * nlem
    ncou = nsph * nsph - 1
    litpo = li + li + 1
    litpos = litpo * litpo
    lmpo = lm + 1
    lmtpo = li + le + 1
    lmtpos = lmtpo * lmtpo
    nv3j = (lm * (lm + 1) * (2 * lm + 7)) / 6;
    ndi = nsph * nlim
    ndit = 2 * ndi
    nllt = 2 * nsph * li * (li + 2) if nlemt == 0 else nlemt
    npnt = gconf['npnt']
    npntts = gconf['npntts']
    max_n = npnt if npnt > npntts else npntts
    nhspo = 2 * max_n - 1
    configurations = sconf['configurations']
    num_layers = 0
    max_layers = 1
    for nli in range(configurations):
        nl = sconf['nshl'][nli]
        if (nli == 0 and gconf['application'] == "INCLU"): nl += 1
        num_layers += nl
        if (nl > max_layers): max_layers = nl
    try:
        max_gpu_ram = int(model['system_settings']['max_gpu_ram'])
        matrix_dim = 2 * gconf['nsph'] * gconf['li'] * (gconf['li'] + 2)
        matrix_size_bytes = 16 * matrix_dim * matrix_dim
        matrix_size_Gb = float(matrix_size_bytes) / 1024.0 / 1024.0 / 1024.0
        print("INFO: estimated matrix size is {0:.3g} GiB.".format(matrix_size_Gb))
        if (max_gpu_ram > 0):
            max_gpu_ram_bytes = max_gpu_ram * 1024 * 1024 * 1024
            if (matrix_size_bytes < max_gpu_ram_bytes):
                max_gpu_processes = int(max_gpu_ram_bytes / matrix_size_bytes)
                print("INFO: system supports up to %d simultaneous processes on GPU."%max_gpu_processes)
                print("INFO: only %d GPU processes allowed, if using refinement."%(max_gpu_processes / 3))
            else:
                print("WARNING: estimated matrix size is larger than available GPU memory!")
        else:
            print("INFO: no GPU RAM declared.")
        max_host_ram = int(model['system_settings']['max_host_ram'])
        max_host_ram_bytes = max_host_ram * 1024 * 1024 * 1024
        if (max_host_ram > 0):
            required_ram_bytes = 0
            if (gconf['application'] == "CLUSTER"):
                # ClusterIterationData section
                required_ram_bytes += 8 * 37
                required_ram_bytes += 8 * 10
                required_ram_bytes += 16
                required_ram_bytes += 4 * 7
                required_ram_bytes += 1
                required_ram_bytes += 8 * (44 + 6 * lm + ndit)
                required_ram_bytes += 8 * (132 + 5 * nsph + 12 * lm)
                required_ram_bytes += 16 * (24 + 4 * nsph + ndit * ndit)
                # ParticleDescriptorCluster section
                required_ram_bytes += 2 * 2
                required_ram_bytes += 4 * 16
                required_ram_bytes += 8 * 20
                required_ram_bytes += 8
                required_ram_bytes += 16 * (2 * li * nsph + 16 + 2 * nhspo + nsph + max_layers + 1)
                required_ram_bytes += 16 * (2 * li + configurations)
                required_ram_bytes += 8 * (num_layers + configurations + 4 * nsph)
                required_ram_bytes += 4 * (nsph + configurations)
                required_ram_bytes += 8 * 12
                required_ram_bytes += 16 * 22 * nsph
                required_ram_bytes += 8 * 7 * nsph
                required_ram_bytes += 8 * 4 * nsph
                required_ram_bytes += 8 * 9
                required_ram_bytes += 16
                required_ram_bytes += 8 * 3
                required_ram_bytes += 16 * (20 + 2 * ndi * nlem + ndit * nlemt)
                required_ram_bytes += 8 * (2 + 2 * ndit)
                required_ram_bytes += 4 * 29
                required_ram_bytes += 8 * 31
                required_ram_bytes += 16
                required_ram_bytes += 16 * (
                    4 * nllt + nlemt * nlemt + 36 + ncou * litpo + nsph * lmtpo
                    + ncou * litpos + nsph * lmtpos
                )
                required_ram_bytes += 4 * ((lm + 1) * lm)
                required_ram_bytes += 8 * (8 + nv3j + lmtpo)
            else:
                print("ERROR: unrecognized application name \"%s\""%gconf['application'])
                raise KeyError("unrecognized application name \"{0:s}\"".format(gconf['application']))
            required_ram_gb = required_ram_bytes / 1024.0 / 1024.0 / 1024.0
            print("INFO: model requires %.5gGiB of host RAM."%required_ram_gb)
            if (required_ram_bytes < max_host_ram_bytes):
                max_host_processes = int(max_host_ram_bytes / required_ram_bytes)
                print("INFO: system supports up to %d simultaneous processes"%max_host_processes)
                print("      (N.B.: not including overheads!)")
            else:
                print("WARNING: estimated matrix size is larger than available host memory!")
        else:
            print("WARNING: no host RAM declared!")
    except KeyError as ex:
        print(ex)
        print("WARNING: missing system description! Cannot estimate recommended execution.")
    cpu_count = multiprocessing.cpu_count()
    print("INFO: the number of detected CPUs is %d."%cpu_count)
    
## \brief Write the geometry configuration dictionary to legacy format.
#
#  \param conf: `dict` Geometry configuration dictionary.
#  \return result: `int` An exit code (0 if successful).
def write_legacy_gconf(conf):
    result = 0
    out_file = str(conf['out_file'])
    nsph = conf['nsph']
    str_line = "INIT"
    # Write legacy output
    output = open(out_file, 'w')
    if (conf['application'] == "SPHERE"):
        str_line = " {0:4d} {1:4d} {2:4d} {3:4d} {4:4d} {5:4d}\n".format(
            nsph, conf['li'], conf['inpol'], conf['npnt'], conf['npntts'], conf['isam']
        )
        output.write(str_line)
    else:
        mxndm = 2 * nsph * conf['li'] * (conf['li'] + 2)
        str_line = " {0:4d} {1:4d} {2:4d} {3:4d} {4:4d} {5:4d} {6:4d} {7:4d} {8:4d}\n".format(
            nsph, conf['li'], conf['le'], mxndm, conf['inpol'],
            conf['npnt'], conf['npntts'], conf['iavm'], conf['isam']
        )
        output.write(str_line)
        for si in range(nsph):
            str_line = " {0:15.8E} {1:15.8E} {2:15.8E}\n".format(
                conf['vec_sph_x'][si], conf['vec_sph_y'][si], conf['vec_sph_z'][si]
            )
            output.write(str_line)
    str_line = " {0:7.2E}  {1:7.2E}  {2:7.2E}  {3:7.2E}  {4:7.2E}  {5:7.2E}\n".format(
        conf['th'], conf['thstp'], conf['thlst'],
        conf['ph'], conf['phstp'], conf['phlst']
    )
    output.write(str_line)
    str_line = " {0:7.2E}  {1:7.2E}  {2:7.2E}  {3:7.2E}  {4:7.2E}  {5:7.2E}\n".format(
        conf['ths'], conf['thsstp'], conf['thslst'],
        conf['phs'], conf['phsstp'], conf['phslst']
    )
    output.write(str_line)
    str_line = " {0:d}\n0\n".format(conf['jwtm'])
    output.write(str_line)
    if (conf['use_refinement']):
        str_line = "USE_REFINEMENT=1\n"
    else:
        str_line = "USE_REFINEMENT=0\n"
    output.write(str_line)
    if (conf['dyn_orders']):
        str_line = "USE_DYN_ORDERS=1\n"
    else:
        str_line = "USE_DYN_ORDERS=0\n"
    output.write(str_line)
    max_host_ram_gb = conf['max_host_ram_gb']
    if (max_host_ram_gb > 0.0):
        str_line = "HOST_RAM_GB={0:.3f}\n".format(max_host_ram_gb)
        output.write(str_line)
    max_gpu_ram_gb = conf['max_gpu_ram_gb']
    if (max_gpu_ram_gb > 0.0):
        str_line = "GPU_RAM_GB={0:.3f}\n".format(max_gpu_ram_gb)
        output.write(str_line)
    output.close()
    return result

## \brief Write the scatterer configuration dictionary to legacy format.
#
#  \param conf: `dict` Scatterer configuration dictionary.
#  \return result: `int` An exit code (0 if successful).
def write_legacy_sconf(conf):
    result = 0
    out_file = str(conf['out_file'])
    nsph = conf['nsph']
    ies = conf['ies']
    exdc = conf['exdc']
    wp = conf['wp']
    xip = conf['xip']
    idfc = conf['idfc']
    instpc = conf['instpc']
    xi_flag = 3
    nxi = conf['nxi']
    # Write legacy output
    output = open(out_file, 'w')
    str_line = " {0:3d}{1:3d}\n".format(nsph, ies)
    output.write(str_line)
    str_line = " {0:12.7E} {1:12.7E} {2:12.7E} {3:2d} {4:4d} {5:4d} {6:3d}\n".format(
        exdc, wp, xip, idfc, nxi, instpc, xi_flag
    )
    output.write(str_line)
    if (instpc == 0):
        for ixi in range(nxi):
            str_line = "{0:.3E}\n".format(conf['vec_xi'][ixi])
            output.write(str_line)
    else:
        str_line = "{0:.3E}  {1:.3E}\n".format(conf['xi_start'], conf['xi_step'])
        output.write(str_line)
    sphere_line_count = 0
    placed_spheres = 0
    last_type = 0
    dedfb_type = 0
    for si in range(nsph):
        if (conf['vec_types'][si] > last_type):
            dedfb_type = placed_spheres + 1
            last_type = conf['vec_types'][si]
        str_line = "{0:5d}".format(dedfb_type)
        output.write(str_line)
        sphere_line_count += 1
        placed_spheres += 1
        if (sphere_line_count == 16):
            output.write("\n")
            sphere_line_count = 0
    if (sphere_line_count != 0):
        output.write("\n")
    for ci in range(conf['configurations']):
        layers = conf['nshl'][ci]
        str_line = "{0:3d}   {1:15.7E}\n".format(layers, conf['ros'][ci])
        output.write(str_line)
        if (conf['application'] == "INCLUSION" and ci == 0):
            layers += 1
        for cj in range(layers):
            str_line = " {0:.7E}\n".format(conf['rcf'][ci][cj])
            output.write(str_line)
    if (conf['application'] != "INCLUSION"):
        if (conf['idfc'] == 0):
            # Write all layers in each configuration for every wavelength
            for xk in range(conf['nxi']):
                for xi in range(conf['configurations']):
                    for xj in range(1 + int(conf['nshl'][xi] / 2)):
                        rdc0 = conf['rdc0'][xj][xi][xk]
                        idc0 = conf['idc0'][xj][xi][xk]
                        if (rdc0 != 0.0 or idc0 != 0.0):
                            str_line = "({0:11.5E},{1:11.5E})\n".format(rdc0, idc0)
                            output.write(str_line)
        elif (conf['idfc'] == -1):
            # Write reference scale constants for each layer in each configuration
            for xi in range(conf['configurations']):
                for xj in range(1 + int(conf['nshl'][xi] / 2)):
                    rdc0 = conf['rdc0'][xj][xi][0]
                    idc0 = conf['idc0'][xj][xi][0]
                    if (rdc0 != 0.0 or idc0 != 0.0):
                        str_line = "({0:11.5E},{1:11.5E})\n".format(rdc0, idc0)
                        output.write(str_line)
    else: # specialized output for INCLUSION
        if (conf['idfc'] == 0):
            # Write all layers in each configuration for every wavelength
            for xk in range(conf['nxi']):
                for xi in range(conf['configurations']):
                    layers = int(conf['nshl'][xi])
                    if (xi == 0):
                        layers += 1
                    for xj in range(layers):
                        rdc0 = conf['rdc0'][xj][xi][xk]
                        idc0 = conf['idc0'][xj][xi][xk]
                        if (rdc0 != 0.0 or idc0 != 0.0):
                            str_line = "({0:11.5E},{1:11.5E})\n".format(rdc0, idc0)
                            output.write(str_line)
        elif (conf['idfc'] == -1):
            # Write reference scale constants for each layer in each configuration
            for xi in range(conf['configurations']):
                layers = int(conf['nshl'][xi])
                if (xi == 0):
                    layers += 1
                for xj in range(layers):
                    rdc0 = conf['rdc0'][xj][xi][0]
                    idc0 = conf['idc0'][xj][xi][0]
                    if (rdc0 != 0.0 or idc0 != 0.0):
                        str_line = "({0:11.5E},{1:11.5E})\n".format(rdc0, idc0)
                        output.write(str_line)
    output.write("0\n")
    output.close()
    return result

## \brief Export the model to a set of OBJ files for 3D visualization.
#
#  This function exports the model as a single OBJ file, containing the
#  information to visualize the particle with 3D software tools. The model
#  file is associated with a MTL material libray file, used to assign colors
#  to spheres of different type.
#
#  \param scatterer: `dict` Scatterer configuration dictionary (gets modified)
#  \param geometry: `dict` Geometry configuration dictionary (gets modified)
#  \param max_rad: `float` Maximum allowed radial extension of the aggregate
def write_obj(scatterer, geometry, max_rad):
    out_dir = scatterer['out_file'].absolute().parent
    out_model_path = Path(str(geometry['out_file']) + ".obj")
    out_material_path = Path(str(geometry['out_file']) + ".mtl")
    # The following color scale is obtained from the normalized
    # RGB values of colors divided by 2.
    color_strings = [
        "0.50 0.50 0.50\n", # white
        "0.30 0.00 0.00\n", # red
        "0.00 0.00 0.30\n", # blue
        "0.00 0.30 0.00\n", # green
        "0.35 0.35 0.00\n", # yellow
        "0.21 0.03 0.34\n", # purple
        "0.00 0.49 0.49\n", # cyan
        "0.27 0.16 0.08\n", # brown
        "0.19 0.21 0.09\n", # olive
        "0.50 0.00 0.50\n", # magenta
        "0.16 0.38 0.37\n", # turquoise
        "0.25 0.25 0.25\n", # gray
        "0.50 0.36 0.38\n", # pink
        "0.50 0.32 0.00\n", # orange
        "0.48 0.45 0.34\n", # vanilla
        "0.16 0.39 0.24\n", # emerald
    ]
    color_names = [
        "white", "red", "blue", "green", "yellow", "purple", "cyan", "brown",
        "olive", "magenta", "turquoise", "gray", "pink", "orange", "vanilla",
        "emerald"
    ]
    mtl_file = open(str(out_material_path), "w")
    for mi in range(len(color_strings)):
        mtl_line = "newmtl "  + color_names[mi] + "\n"
        mtl_file.write(mtl_line)
        color_line = color_strings[mi]
        mtl_file.write("   Ka " + color_line)
        mtl_file.write("   Ks " + color_line)
        mtl_file.write("   Kd " + color_line)
        mtl_file.write("   Ns 100.0\n")
        mtl_file.write("   Tr 0.0\n")
        mtl_file.write("   illum 2\n\n")
    mtl_file.close()
    pl = pv.Plotter()
    for si in range(scatterer['nsph']):
        sph_type_index = scatterer['vec_types'][si]
        # color_index = 1 + (sph_type_index % (len(color_strings) - 1))
        # color_by_name = color_names[sph_type_index]
        radius = scatterer['ros'][sph_type_index - 1] / max_rad
        x = geometry['vec_sph_x'][si] / max_rad
        y = geometry['vec_sph_y'][si] / max_rad
        z = geometry['vec_sph_z'][si] / max_rad
        mesh = pv.Sphere(radius, (x, y, z))
        pl.add_mesh(mesh, color=None)
    pl.export_obj(str(Path(str(out_dir), "TMP_MODEL.obj")))
    tmp_model_file = open(str(Path(str(out_dir), "TMP_MODEL.obj")), "r")
    out_model_file = open(str(out_model_path), "w")
    mtl_line = "mtllib {0:s}\n".format(out_material_path.name)
    sph_index = 0
    sph_type_index = 0
    old_sph_type_index = 0
    str_line = tmp_model_file.readline()
    while (str_line != ""):
        if (str_line.startswith("mtllib")):
            str_line = mtl_line
        elif (str_line.startswith("g ")):
            sph_index += 1
            sph_type_index = scatterer['vec_types'][sph_index - 1]
            if (sph_type_index == old_sph_type_index):
                str_line = tmp_model_file.readline()
                str_line = tmp_model_file.readline()
            else:
                old_sph_type_index = sph_type_index
                color_index = 1 + (sph_type_index - 1) % (len(color_names) - 1)
                str_line = "g grp{0:04d}\n".format(sph_type_index - 1)
                out_model_file.write(str_line)
                str_line = tmp_model_file.readline()
                str_line = "usemtl {0:s}\n".format(color_names[color_index])
        out_model_file.write(str_line)
        str_line = tmp_model_file.readline()
    out_model_file.close()
    tmp_model_file.close()
    os.remove(str(Path(str(out_dir), "TMP_MODEL.obj")))
    os.remove(str(Path(str(out_dir), "TMP_MODEL.mtl")))

## \brief Exit code (0 for success)
exit_code = main()
exit(exit_code)
