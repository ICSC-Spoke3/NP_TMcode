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

## @package inertia
#  \brief Compute the momentum of inertia for a model particle.
#
#  The complete description of a model particle requires knowledge of some properties,
#  such as filling factor and symmetry. These can be estimated by comparing the
#  model particle with its ellipsoiod of inertia. The purpose of this script is to
#  perform such comparison on NP_TMcode model particles.

import math
import numpy as np
#import pdb

from sys import argv

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
        print("\nType \"./inertia.py --help\" to get more detailed help.")
        errors = 1
    if config['help_mode']:
        config['help_mode'] = True
        print_help()
    else:
        if config['scat_name'] == '':
            print("ERROR: no scatterer configuration file given (missing --scat).")
            errors += 1
        if config['geom_name'] == '':
            print("ERROR: no scatterer configuration file given (missing --geom).")
            errors += 1
        if errors == 0:
            sph_types = get_types(config)
            centers = get_centers(config)
            if (sph_types['n'] != centers['n']):
                raise Exception("ERROR: Geometry file has not the same number of spheres as scatterer file!")
            # INFO CODE SECTION
            num_types = len(sph_types['num_layers'])
            for i in range(num_types):
                type_mass = get_sphere_type_mass(i + 1, sph_types, config)
                print("INFO: mass of type %d is %.5e kg."%(i + 1, type_mass))
            masses = np.array([
                get_sphere_mass(i + 1, sph_types, config) for i in range(len(sph_types['vec_types']))
            ])
            particle_mass = get_particle_mass(sph_types, config)
            print("INFO: total particle mass is %.5e kg."%(particle_mass))
            Itot = get_inertia_tensor(masses, sph_types, centers)
            print("INFO: Itot =", Itot)
            Itot_cm, cm = get_cm_inertia_tensor(masses, sph_types, centers)
            print("INFO: cm =", cm)
            print("INFO: Itot_cm =", Itot_cm)
            fill_factor = filling_factor(sph_types, Itot_cm, config)
            print("INFO: the particle's filling factor is %.5g"%fill_factor)
            # END OF INFO CODE SECTION
    return errors

## \brief Compute the particle's filling factor.
#
#  In this application the filling factor is defined as the ratio of the
#  sum of the volumes of the spheres that compose the particle with respect
#  to the volume of a homogeneous ellipsoid that has the same rotational
#  properties as the particle.
#
#  \param stypes: `dict` Dictionary of sphere types.
#  \param i_tot_cm: `numpy.ndarray` Inertia tensor of the particle with respect
#  to the center of mass.
#  \param config: `dict` Configuration dictionary.
#  \return ff: `float` The particle's filling factor.
def filling_factor(stypes, i_tot_cm, config):
    ff = 0.0
    vtot = get_total_volume(stypes, config)
    mtot = get_particle_mass(stypes, config)
    eigenval, eigenvec = np.linalg.eigh(i_tot_cm)
    print("INFO: ellipsoid axes:", eigenval)
    print("INFO: ellipsoid directions:", eigenvec)
    gf = 5.0 / (2.0 * mtot)
    I1, I2, I3 = eigenval
    a = np.sqrt(max(0.0, gf * (I2 + I3 - I1)))
    b = np.sqrt(max(0.0, gf * (I1 + I3 - I2)))
    c = np.sqrt(max(0.0, gf * (I1 + I2 - I3)))
    # Volume of the ellipsoid that has the same moments of inertia as the particle
    vell = 4.0 * math.pi * a * b * c / 3.0
    ff = vtot / vell
    return ff

## \brief Create a dictionary of sphere types.
#
#  This function opens a NP_TMcode scatterer configuration file and creates
#  sphere type dictionary from it.
#
#  \param config: `dict` Configuration dictionary.
#  \return centers: `dict` Sphere centers dictionary.
def get_centers(config):
    centers = {
        'n': 0,
        'x': [],
        'y': [],
        'z': []
    }
    geom_file = open(config['geom_name'], 'r')
    fline = geom_file.readline()
    elems = ingest_line(fline)
    nsph = int(elems[0])
    centers['n'] = nsph
    for i in range(nsph):
        fline = geom_file.readline()
        elems = ingest_line(fline)
        value = float(elems[0])
        centers['x'].append(value)
        value = float(elems[1])
        centers['y'].append(value)
        value = float(elems[2])
        centers['z'].append(value)
    geom_file.close()
    return centers

## \brief Compute the tensor of inertia of an aggregate with respect to the center of mass.
#
#  \param masses: `numpy.ndarray` Array of masses in the aggregate.
#  \param stypes: `dict` Dictionary of sphere types.
#  \param centers: `dict` Dictionary containing the positions of the spheres.
#  \return Itot: `numpy.ndarray` The total tensor of inertia with respect to the center of mass.
def get_cm_inertia_tensor(masses, stypes, centers):
    N = len(masses)
    Itot = np.zeros((3, 3))
    
    m_tot = np.sum(masses)
    positions = np.array([
        [centers['x'][i], centers['y'][i], centers['z'][i]] for i in range(len(centers['x']))
    ])
    center_of_mass = np.sum(masses[:, np.newaxis] * positions, axis=0) / m_tot
    for i in range(N):
        M = masses[i]
        typeID = stypes['vec_types'][i] - 1
        R = stypes['radii'][typeID]
        x = centers['x'][i] - center_of_mass[0]
        y = centers['y'][i] - center_of_mass[1]
        z = centers['z'][i] - center_of_mass[2]
        
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
    return Itot, center_of_mass

## \brief Compute the tensor of inertia of an aggregate with respect to origin of coordinates.
#
#  \param masses: `list-like` Array of masses in the aggregate.
#  \param stypes: `dict` Dictionary of sphere types.
#  \param centers: `dict` Dictionary containing the positions of the spheres.
#  \return Itot: `numpy.ndarray` The total tensor of inertia with respect to the reference frame.
def get_inertia_tensor(masses, stypes, centers):
    N = len(masses)
    Itot = np.zeros((3, 3))

    for i in range(N):
        M = masses[i]
        typeID = stypes['vec_types'][i] - 1
        R = stypes['radii'][typeID]
        x = centers['x'][i]
        y = centers['y'][i]
        z = centers['z'][i]
        
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
    return Itot

## \brief Compute the total particle mass.
#
#  \param stypes: `dict` Sphere types dictionary.
#  \param config: `dict` Configuration dictionary.
#  \return mass: `float` The mass of the particle, expressed in kg.
def get_particle_mass(stypes, config):
    # itype = stypes['vec_types'][isphere - 1]
    mass = 0.0
    for i in range(len(stypes['vec_types'])):
        mass += get_sphere_mass(i + 1, stypes, config)
    return mass

## \brief Compute the mass of a given sphere.
#
#  \param isphere: `int` ID of the sphere in the aggregate (base 1).
#  \param stypes: `dict` Sphere types dictionary.
#  \param config: `dict` Configuration dictionary.
#  \return mass: `float` The mass of the sphere, expressed in kg.
def get_sphere_mass(isphere, stypes, config):
    itype = stypes['vec_types'][isphere - 1]
    mass = get_sphere_type_mass(itype, stypes, config)
    return mass

## \brief Compute the mass of a given type of sphere.
#
#  \param itype: `int` ID of the sphere type (base 1).
#  \param stypes: `dict` Sphere types dictionary.
#  \param config: `dict` Configuration dictionary.
#  \return type_mass: `float` The mass of the type of sphere, expressed in kg.
def get_sphere_type_mass(itype, stypes, config):
    type_mass = 0.0
    materials_vec = stypes['material_ids'][itype - 1]
    radius = stypes['radii'][itype - 1]
    same_material = True
    for i in range(1, len(materials_vec)):
        if (materials_vec[i] != materials_vec[0]):
            same_material = False
            break # for i
    # end of i loop
    if same_material:
        # Compute the mass of a homogeneous sphere
        material_id = materials_vec[0]
        specw = config['specific'][material_id - 1] # this value is in g / cm^3
        volume = 4.0 * math.pi / 3.0 * math.pow(radius, 3.0)
        type_mass = volume * specw * 1.0e3 # this value is in kg
    else:
        # Compute the mass of a multi-layered sphere
        num_layers = stypes['num_layers'][itype - 1]
        ros = stypes['radii'][itypes - 1]
        rcf = stypes['frac_radii'][itypes - 1]
        radii = [ros * rcf[i] for i in range(num_layers)]
        inner_radius = 0.0
        for i in range(num_layers):
            layer_index = 1 + (i + 1) / 2
            specw = 0.0
            if (layer_index % 2 == 1):
                material_id = materials_vec[layer_index - 1]
                specw = config['specific'][material_id - 1] # this value is in g / cm^3
            else:
                # compute average specw
                material_1 = materials_vec[i - 1]
                specw_1 = config['specific'][material_1 - 1] # this value is in g / cm^3
                material_2 = materials_vec[i]
                specw_2 = config['specific'][material_2 - 1] # this value is in g / cm^3
                specw = (spwcw_1 + specw_2 ) / 2.0
            volume = math.pow(radii[i], 3.0)
            volume -= math.pow(inner_radius, 3.0)
            volume *= (4.0 * math.pi / 3.0)
            inner_radius = radii[i]
            type_mass += volume * specw * 1.0e3 # this value is in kg
            # i loop ends here
        # Compute the mass of coating, if present
        if (itype == 1 and config['ies'] != 0):
            coat_id = len(stypes['frac_radii'][0]) - 1
            coat_radius = stypes['radii'][0] * stypes['frac_radii'][0][coat_id]
            coat_volume = math.pow(coat_radius, 3.0)
            for i in range(1, stypes['vec_types']):
                radius = stypes['radii'][stypes['vec_types'][i] - 1]
                coat_volume -= math.pow(radius, 3.0)
                # i loop ends here
            if (coat_volume < 0.0):
                raise Exception("ERROR: negative coating volume!")
            coat_volume *= (4.0 * math.pi / 3.0)
            coat_specw_id = stypes['material_ids'][0][len(stypes['material_ids'][0]) - 1]
            coat_specw = config['specific'][coat_specw_id - 1] # this value is in g / cm^3
            type_mass += coat_volume * coat_specw * 1.0e3 # this value is in kg 
    return type_mass

## \brief Compute the total volume of the spheres.
#
#  \param stypes: `dict` Sphere types dictionary.
#  \param config: `dict` Configuration dictionary.
#  \return vtot: `float` The sum of the volumes of the spheres.
def get_total_volume(stypes, config):
    vtot = 0.0
    if config['ies'] == 0:
        for i in range(len(stypes['vec_types'])):
            typeID = stypes['vec_types'][i] - 1
            radius = stypes['radii'][typeID]
            vtot += math.pow(radius, 3.0)
    else:
        factor = stypes['frac_radii'][0][len(stypes['frac_radii'][0]) - 1]
        radius = factor * stypes['radii'][0]
        vtot = math.pow(radius, 3.0)
    return (4.0 * math.pi / 3.0) * vtot

## \brief Create a dictionary of sphere types.
#
#  This function opens a NP_TMcode scatterer configuration file and creates
#  sphere type dictionary from it.
#
#  \param config: `dict` Configuration dictionary.
#  \return stypes: `dict` Sphere types dictionary.
def get_types(config):
    stypes = {
        'n': 0,
        'num_types': 0,
        'vec_types': [],
        'num_layers': [],
        'radii': [],
        'frac_radii': [],
        'material_ids': []
    }
    sc_file = open(config['scat_name'], 'r')
    fline = sc_file.readline() # Get NSPH and IES
    split_line = ingest_line(fline)
    config['nsph'] = int(split_line[0])
    config['ies'] = int(split_line[1])
    stypes['n'] = config['nsph']
    fline = sc_file.readline() # Get NXI
    split_line = ingest_line(fline)
    nxi = int(split_line[4])
    instpc = int(split_line[5])
    if (instpc == 1):
        nxi = 2
    for li in range(nxi):
        fline = sc_file.readline()
    # Detect the number of types
    read_spheres = 0
    type_id = 0
    while (read_spheres < config['nsph']):
        fline = sc_file.readline()
        split_line = ingest_line(fline)
        for ti in range(len(split_line)):
            cur_type = int(split_line[ti])
            if (cur_type > read_spheres + ti):
                type_id += 1
            if (len(stypes['vec_types']) < config['nsph']):
                stypes['vec_types'].append(type_id)
        read_spheres += len(split_line)
    stypes['num_types'] = type_id
    # Detect the layers and the radii of each type
    read_types = 0
    while (read_types < stypes['num_types']):
        # read a type
        fline = sc_file.readline()
        split_line = ingest_line(fline)
        num_layers = int(split_line[0])
        radius = float(split_line[1].replace('d', 'e').replace('D', 'E'))
        stypes['num_layers'].append(num_layers)
        stypes['radii'].append(radius)
        stypes['frac_radii'].append([])
        stypes['material_ids'].append([])
        if (read_types == 0 and config['ies'] != 0):
            num_layers += 1
        for li in range(num_layers):
            fline = sc_file.readline()
            frac_rad = float(fline.replace('d', 'e').replace('D', 'E'))
            stypes['frac_radii'][read_types].append(frac_rad)
        read_types += 1
    # Detect the material IDs
    last_dc0 = (0.0,0.0)
    read_types = 0
    found_dc0s = []
    num_dc0s = 0
    while (read_types < stypes['num_types']):
        read_layers = 0
        material_id = 0
        num_layers = int(stypes['num_layers'][read_types] / 2) + 1
        if (read_types == 0 and config['ies'] != 0):
            num_layers += 1
        while (read_layers < num_layers):
            fline = sc_file.readline()
            split_line = fline.split('(')[1].split(',')
            rval = float(split_line[0])
            ival = float(split_line[1][:-2])
            dc0 = (rval, ival)
            if not dc0 in found_dc0s:
                found_dc0s.append(dc0)
                num_dc0s += 1
                stypes['material_ids'][read_types].append(num_dc0s)
            else:
                material_id = 1 + index_of(dc0, found_dc0s)
                stypes['material_ids'][read_types].append(material_id)
            read_layers += 1
        read_types += 1
    print("DEBUG: stypes =", stypes)
    sc_file.close()
    return stypes

## \brief Get the index of an element in an array.
#
#  \param element: Anything.
#  \param array: `array-like` Array with the same type of the element sought for.
#  \return index: `int` 0-based index of the element.
def index_of(element, array):
    index = 0
    while (element != array[index]):
        index += 1
    return index
           
## \brief Tranform a line in a sequence of string elements.
#
#  \param line: `string` Input file line
#  \return result: List of `string` elements.
def ingest_line(line):
    result = []
    elem = ""
    for c in line:
        if c != ' ' and c != '\t' and c != '\n' :
            elem += c
        else:
            if elem != "":
                result.append(elem)
                elem = ""
    if elem != "":
        result.append(elem)
    return result

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
        'scat_name': '',
        'geom_name': '',
        'specific': [1.0],
        'help_mode': False
    }
    arg_index = 1
    skip_arg = False
    for arg in argv[1:]:
        if skip_arg:
            skip_arg = False
            continue
        split_arg = arg.split('=')
        if (arg.startswith("--scat")):
            config['scat_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--geom")):
            config['geom_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--spew=")):
            #config['specific'] = float(split_arg[1])
            config['specific'] = []
            str_specifics = split_arg[1]
            vec_specifics = str_specifics.split(',')
            for si in range(len(vec_specifics)):
                config['specific'].append(float(vec_specifics[si]))
        elif (arg.startswith("--help")):
            config['help_mode'] = True
        else:
            raise ValueError("Unrecognized argument \'{0:s}\'".format(arg))
        arg_index += 1
    return config

## \brief Print a command-line help summary.
def print_help():
    print("                                                                ")
    print("***                COMPUTE MOMENTUM OF INERTIA               ***")
    print("                                                                ")
    print("Compute the momentum of inertia for a NP_TMcode particle model. ")
    print("                                                                           ")
    print("Usage: \"./inertia.py --scat SCAT --geom GEOM [OPTIONS]\"                  ")
    print("                                                                           ")
    print("Valid options are:                                                         ")
    print("--scat SCAT              Scatterer configuration file (mandatory).         ")
    print("--geom GEOM              Geometry configuration file (mandatory).          ")
    print("--help                   Print this help and exit.                         ")
    print("--spew=SPEC_WEIGHT       Specific weight of materials (in g/cm3, optional).")
    print("                                                                           ")

# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
exit(res)
