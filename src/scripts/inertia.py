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
#  \ingroup python_tools
#  \brief Compute the mass, the tensor of inertia and the filling factor of a model particle.
#
#  The complete description of a model particle requires knowledge of some
#  properties, such as filling factor and symmetry. These can be estimated by
#  comparing the model particle with its ellipsoid of inertia. The script
#  reads the scatterer (`DEDFB`) and geometry (`DSPH`, `DCLU` or `DINCLU`)
#  files of a model and computes:
#
#  - the mass of each sphere type and of the particle;
#  - the tensor of inertia with respect to the origin of the coordinates and
#    to the centre of mass;
#  - the semi-axes of the homogeneous ellipsoid with the same mass and
#    principal moments of inertia, and the filling factor (volume of the
#    particle over the volume of that ellipsoid).
#
#  The model is the one of the code:
#
#  - IOG(i) is the index of the first sphere with the same type as sphere i;
#  - a sphere type has NSHL layers with outer radii RCF(j) * ROS. The odd
#    layers are made of the materials declared for the type, in order; the
#    even layers are transition layers, whose composition changes from the
#    material inside to the one outside with the profile of the code
#    (3 x^2 - 2 x^3, x going from 0 to 1 across the layer). The material
#    outside the last layer is the external medium when NSHL is even;
#  - with IES = 1 (`np_inclusion`) the first type has one more radius and one
#    more material: the external sphere, centred in the origin of the
#    coordinates, made of that material and containing all the spheres.
#
#  The materials are identified by their dielectric constants at the first
#  wavelength of the scatterer file, in order of first appearance; `--spew`
#  gives their specific weights in the same order.

import math
import re
import numpy as np
from numpy.polynomial import Polynomial

from sys import argv

## \cond
from nptm_version import __version__
## \endcond

## \brief Main execution code.
#
# `main()` is the function that handles the creation of the script configuration
# and the execution of the calculation. It returns an integer value corresponding
# to the script exit code.
#
# \returns errors: `int` Number of runtime errors (0 means successful run).
def main():
    errors = 0
    try:
        config = parse_arguments()
    except ValueError as ex:
        print("ERROR: %s"%str(ex))
        print("\nType \"./inertia.py --help\" to get more detailed help.")
        return 1
    if config['version_mode']:
        print("inertia.py v%s."%__version__)
        return 0
    if config['help_mode']:
        print_help()
        return 0
    if config['scat_name'] == '':
        print("ERROR: no scatterer configuration file given (missing --scat).")
        errors += 1
    if config['geom_name'] == '':
        print("ERROR: no geometry configuration file given (missing --geom).")
        errors += 1
    if errors > 0:
        return errors
    try:
        stypes = get_types(config)
        centers = get_centers(config)
        if (stypes['n'] != centers['n']):
            raise ValueError(
                "the geometry file has %d spheres, the scatterer file %d."%(centers['n'], stypes['n'])
            )
        specw = get_specific_weights(stypes, config)
    except (OSError, ValueError, IndexError) as ex:
        print("ERROR: %s"%str(ex))
        return 1
    print("INFO: %d spheres of %d types, %d materials%s."%(
        stypes['n'], stypes['num_types'], stypes['num_materials'],
        ", with an external sphere" if stypes['ies'] != 0 else ""
    ))
    for mi in range(stypes['num_materials']):
        dc0 = stypes['materials'][mi]
        print("INFO: material %d (dielectric constant (%.5g, %.5g)) has specific weight %.5g g/cm^3."%(
            mi + 1, dc0.real, dc0.imag, specw[mi]
        ))
    for ti in range(stypes['num_types']):
        mass, moment = get_type_mass_moment(ti + 1, stypes, specw, config)
        print("INFO: mass of type %d is %.5e kg."%(ti + 1, mass))
    bodies = get_bodies(stypes, centers, specw, config)
    particle_mass = sum(b[0] for b in bodies)
    print("INFO: total particle mass is %.5e kg."%(particle_mass))
    if (particle_mass <= 0.0):
        print("ERROR: the particle mass is not positive: check the specific weights.")
        return 1
    np.set_printoptions(precision=5)
    Itot = get_inertia_tensor(bodies)
    print("INFO: tensor of inertia with respect to the origin (kg m^2):\n", Itot)
    Itot_cm, cm = get_cm_inertia_tensor(bodies)
    print("INFO: centre of mass (m):", cm)
    print("INFO: tensor of inertia with respect to the centre of mass (kg m^2):\n", Itot_cm)
    fill_factor, axes, directions = filling_factor(stypes, centers, Itot_cm, particle_mass)
    print("INFO: semi-axes of the equivalent ellipsoid (m):", axes)
    print("INFO: directions of the semi-axes (columns):\n", directions)
    print("INFO: the particle's filling factor is %.5g"%fill_factor)
    return 0

## \brief Compute the particle's filling factor.
#
#  In this application the filling factor is defined as the ratio of the
#  volume of the particle (the sum of the volumes of the spheres, or the
#  volume of the external sphere for an inclusion) to the volume of the
#  homogeneous ellipsoid that has the same mass and the same principal
#  moments of inertia as the particle. For an ellipsoid of mass M and
#  semi-axes a, b, c, the principal moments are I1 = M (b^2 + c^2) / 5 and
#  the like, so a^2 = 5 (I2 + I3 - I1) / (2 M).
#
#  \param stypes: `dict` Dictionary of sphere types.
#  \param centers: `dict` Dictionary of sphere centres.
#  \param i_tot_cm: `numpy.ndarray` Inertia tensor of the particle with respect
#  to the centre of mass.
#  \param mtot: `float` Mass of the particle.
#  \return result: `tuple(float, numpy.ndarray, numpy.ndarray)` The filling
#  factor, the semi-axes of the ellipsoid and their directions (columns).
def filling_factor(stypes, centers, i_tot_cm, mtot):
    vtot = get_total_volume(stypes)
    eigenval, eigenvec = np.linalg.eigh(i_tot_cm)
    gf = 5.0 / (2.0 * mtot)
    I1, I2, I3 = eigenval
    a = math.sqrt(max(0.0, gf * (I2 + I3 - I1)))
    b = math.sqrt(max(0.0, gf * (I1 + I3 - I2)))
    c = math.sqrt(max(0.0, gf * (I1 + I2 - I3)))
    vell = 4.0 * math.pi * a * b * c / 3.0
    ff = vtot / vell if vell > 0.0 else float('inf')
    return (ff, np.array([a, b, c]), eigenvec)

## \brief Read the sphere centres from the geometry file.
#
#  A `DSPH` file (6 values on the first line) describes a single sphere in
#  the origin; `DCLU` and `DINCLU` files list the coordinates of the spheres
#  after the first line.
#
#  \param config: `dict` Configuration dictionary.
#  \return centers: `dict` Sphere centres dictionary (`n`, `x`, `y`, `z`, in meters).
def get_centers(config):
    centers = {
        'n': 0,
        'x': [],
        'y': [],
        'z': []
    }
    with open(config['geom_name'], 'r') as geom_file:
        elems = geom_file.readline().split()
        nsph = int(elems[0])
        centers['n'] = nsph
        if (len(elems) <= 6):
            # DSPH: one sphere in the origin
            for i in range(nsph):
                centers['x'].append(0.0)
                centers['y'].append(0.0)
                centers['z'].append(0.0)
        else:
            for i in range(nsph):
                elems = geom_file.readline().split()
                if (len(elems) < 3):
                    raise ValueError("sphere %d has no coordinates in %s."%(i + 1, config['geom_name']))
                centers['x'].append(to_float(elems[0]))
                centers['y'].append(to_float(elems[1]))
                centers['z'].append(to_float(elems[2]))
    return centers

## \brief Convert a FORTRAN-style number to a float.
#
#  \param text: `str` The number, possibly with a `D` exponent.
#  \return value: `float` The number.
def to_float(text):
    return float(text.replace('D', 'E').replace('d', 'E'))

## \brief Read the sphere types from the scatterer file.
#
#  The file is read as the code reads it (`ScattererConfiguration`): the
#  header, the wavelength section, IOG, then for each type NSHL, ROS and the
#  RCF values (one more for the first type with IES = 1), then the
#  dielectric constants of the first wavelength: (NSHL + 1) / 2 per type,
#  one more for the first type with IES = 1.
#
#  \param config: `dict` Configuration dictionary.
#  \return stypes: `dict` Sphere types dictionary: `n` (spheres), `ies`,
#  `num_types`, `vec_types` (type of each sphere, base 1), and for each
#  type `num_layers`, `radii` (ROS), `frac_radii` (RCF) and `material_ids`
#  (base 1, one per declared dielectric constant); `materials` lists the
#  distinct dielectric constants, `num_materials` their number.
def get_types(config):
    stypes = {
        'n': 0,
        'ies': 0,
        'num_types': 0,
        'vec_types': [],
        'num_layers': [],
        'radii': [],
        'frac_radii': [],
        'material_ids': [],
        'materials': [],
        'num_materials': 0
    }
    with open(config['scat_name'], 'r') as sc_file:
        lines = [line.strip() for line in sc_file if line.strip() != '']
    li = 0
    split_line = lines[li].split()
    nsph = int(split_line[0])
    ies = int(split_line[1])
    stypes['n'] = nsph
    stypes['ies'] = ies
    li += 1
    split_line = lines[li].split()
    nxi = int(split_line[4])
    instpc = int(split_line[5])
    li += 1
    # One line with the first scale and the step, or one line per scale.
    li += 1 if instpc != 0 else nxi
    # IOG(i): index of the first sphere of the same type.
    iog = []
    while (len(iog) < nsph):
        iog += [int(value) for value in lines[li].split()]
        li += 1
    type_id = 0
    for i in range(nsph):
        if (iog[i] >= i + 1):
            type_id += 1
            stypes['vec_types'].append(type_id)
        else:
            stypes['vec_types'].append(stypes['vec_types'][iog[i] - 1])
    stypes['num_types'] = type_id
    # Layers and radii of each type.
    for ti in range(stypes['num_types']):
        split_line = lines[li].split()
        li += 1
        num_layers = int(split_line[0])
        stypes['num_layers'].append(num_layers)
        stypes['radii'].append(to_float(split_line[1]))
        nsh = num_layers + (ies if ti == 0 else 0)
        stypes['frac_radii'].append([to_float(lines[li + j].split()[0]) for j in range(nsh)])
        li += nsh
    # Dielectric constants of the first wavelength.
    dc0_reg = re.compile(r'\(\s*([^,\s]+)\s*,\s*([^)\s]+)\s*\)')
    for ti in range(stypes['num_types']):
        ici = (stypes['num_layers'][ti] + 1) // 2 + (ies if ti == 0 else 0)
        stypes['material_ids'].append([])
        for ic in range(ici):
            match = dc0_reg.search(lines[li])
            if (match is None):
                raise ValueError("expected a dielectric constant \"(re,im)\" at \"%s\" in %s."%(
                    lines[li], config['scat_name']
                ))
            li += 1
            dc0 = complex(to_float(match.group(1)), to_float(match.group(2)))
            if not dc0 in stypes['materials']:
                stypes['materials'].append(dc0)
            stypes['material_ids'][ti].append(1 + stypes['materials'].index(dc0))
    stypes['num_materials'] = len(stypes['materials'])
    return stypes

## \brief Specific weight of each material.
#
#  \param stypes: `dict` Sphere types dictionary.
#  \param config: `dict` Configuration dictionary.
#  \return specw: `list(float)` Specific weight of each material, in g/cm^3.
def get_specific_weights(stypes, config):
    values = config['specific']
    if (len(values) == 1):
        return [values[0] for i in range(stypes['num_materials'])]
    if (len(values) != stypes['num_materials']):
        raise ValueError(
            "--spew gives %d specific weights, but the model has %d materials."%(
                len(values), stypes['num_materials'])
        )
    return values

## \brief Mass and moment of inertia of a sphere type.
#
#  The moment of inertia J is the one of the sphere about any axis through
#  its centre (J = 2 M R^2 / 5 for a homogeneous sphere): each layer adds
#  M = 4 pi Int rho r^2 dr and J = 8 pi / 3 Int rho r^4 dr. For the first
#  type with IES = 1, the external sphere is not included (see
#  `get_bodies()`).
#
#  \param itype: `int` ID of the sphere type (base 1).
#  \param stypes: `dict` Sphere types dictionary.
#  \param specw: `list(float)` Specific weights of the materials, in g/cm^3.
#  \param config: `dict` Configuration dictionary.
#  \return result: `tuple(float, float)` Mass (kg) and moment of inertia (kg m^2).
def get_type_mass_moment(itype, stypes, specw, config):
    ti = itype - 1
    ros = stypes['radii'][ti]
    rcf = stypes['frac_radii'][ti]
    material_ids = stypes['material_ids'][ti]
    num_layers = stypes['num_layers'][ti]
    # Densities in kg/m^3 of the declared materials of the type; the one
    # after the last is the external medium, or the material of the external
    # sphere for the first type of an inclusion.
    rho = [1.0e3 * specw[mid - 1] for mid in material_ids]
    rho_out = 1.0e3 * config['medium']
    mass = 0.0
    moment = 0.0
    r_in = 0.0
    for ns in range(1, num_layers + 1):
        r_out = rcf[ns - 1] * ros
        if (ns % 2 == 1):
            ic = (ns + 1) // 2
            layer_mass, layer_moment = shell_mass_moment(r_in, r_out, rho[ic - 1], rho[ic - 1])
        else:
            ic = ns // 2
            rho_next = rho[ic] if ic < len(rho) else rho_out
            layer_mass, layer_moment = shell_mass_moment(r_in, r_out, rho[ic - 1], rho_next)
        mass += layer_mass
        moment += layer_moment
        r_in = r_out
    return (mass, moment)

## \brief Mass and moment of inertia of a spherical shell.
#
#  The density goes from `rho_in` at `r_in` to `rho_out` at `r_out` with the
#  profile of the transition layers of the code, 3 x^2 - 2 x^3; it is
#  constant if the two values are equal.
#
#  \param r_in: `float` Inner radius (m).
#  \param r_out: `float` Outer radius (m).
#  \param rho_in: `float` Density at the inner radius (kg/m^3).
#  \param rho_out: `float` Density at the outer radius (kg/m^3).
#  \return result: `tuple(float, float)` Mass (kg) and moment of inertia about an axis through the centre (kg m^2).
def shell_mass_moment(r_in, r_out, rho_in, rho_out):
    if (r_out <= r_in):
        return (0.0, 0.0)
    # Integrate in x = (r - r_in) / width, from 0 to 1: with r as variable,
    # a thin layer would lose all its digits in the difference of two large
    # values.
    width = r_out - r_in
    x = Polynomial([0.0, 1.0])
    rho = rho_in + (rho_out - rho_in) * (3.0 * x**2 - 2.0 * x**3)
    r = Polynomial([r_in, width])
    mass_integral = (rho * r**2).integ()
    moment_integral = (rho * r**4).integ()
    mass = 4.0 * math.pi * width * (mass_integral(1.0) - mass_integral(0.0))
    moment = 8.0 * math.pi / 3.0 * width * (moment_integral(1.0) - moment_integral(0.0))
    return (mass, moment)

## \brief The bodies that compose the particle.
#
#  Each body is a sphere described by its mass, its centre and its moment of
#  inertia about an axis through the centre. For an inclusion, the external
#  sphere is a full sphere of its material in the origin, plus one body of
#  negative mass for each sphere it contains (the volume that its material
#  does not fill).
#
#  \param stypes: `dict` Sphere types dictionary.
#  \param centers: `dict` Sphere centres dictionary.
#  \param specw: `list(float)` Specific weights of the materials, in g/cm^3.
#  \param config: `dict` Configuration dictionary.
#  \return bodies: `list(tuple(float, numpy.ndarray, float))` Mass, centre and moment of each body.
def get_bodies(stypes, centers, specw, config):
    bodies = []
    type_values = [
        get_type_mass_moment(ti + 1, stypes, specw, config) for ti in range(stypes['num_types'])
    ]
    for i in range(stypes['n']):
        mass, moment = type_values[stypes['vec_types'][i] - 1]
        center = np.array([centers['x'][i], centers['y'][i], centers['z'][i]])
        bodies.append((mass, center, moment))
    if (stypes['ies'] != 0):
        r_ext = stypes['radii'][0] * stypes['frac_radii'][0][-1]
        rho_ext = 1.0e3 * specw[stypes['material_ids'][0][-1] - 1]
        mass, moment = shell_mass_moment(0.0, r_ext, rho_ext, rho_ext)
        bodies.append((mass, np.zeros(3), moment))
        for i in range(stypes['n']):
            ti = stypes['vec_types'][i] - 1
            radius = stypes['radii'][ti] * stypes['frac_radii'][ti][stypes['num_layers'][ti] - 1]
            mass, moment = shell_mass_moment(0.0, radius, rho_ext, rho_ext)
            center = np.array([centers['x'][i], centers['y'][i], centers['z'][i]])
            bodies.append((-mass, center, -moment))
    return bodies

## \brief Compute the tensor of inertia of the particle with respect to a point.
#
#  Each body adds its own moment on the diagonal and, by the parallel axis
#  theorem, M (|d|^2 E - d d^T), with d the position of its centre with
#  respect to the point.
#
#  \param bodies: `list` Bodies of the particle (see `get_bodies()`).
#  \param origin: `numpy.ndarray` The point (optional, default: origin of the coordinates).
#  \return Itot: `numpy.ndarray` The tensor of inertia (kg m^2).
def get_inertia_tensor(bodies, origin=None):
    if (origin is None):
        origin = np.zeros(3)
    Itot = np.zeros((3, 3))
    for mass, center, moment in bodies:
        d = center - origin
        Itot += moment * np.eye(3)
        Itot += mass * (np.dot(d, d) * np.eye(3) - np.outer(d, d))
    return Itot

## \brief Compute the tensor of inertia of the particle with respect to its centre of mass.
#
#  \param bodies: `list` Bodies of the particle (see `get_bodies()`).
#  \return result: `tuple(numpy.ndarray, numpy.ndarray)` The tensor of
#  inertia (kg m^2) and the centre of mass (m).
def get_cm_inertia_tensor(bodies):
    m_tot = sum(b[0] for b in bodies)
    center_of_mass = sum(b[0] * b[1] for b in bodies) / m_tot
    return (get_inertia_tensor(bodies, center_of_mass), center_of_mass)

## \brief Compute the volume of the particle.
#
#  The sum of the volumes of the spheres or, for an inclusion, the volume of
#  the external sphere.
#
#  \param stypes: `dict` Sphere types dictionary.
#  \return vtot: `float` The volume (m^3).
def get_total_volume(stypes):
    vtot = 0.0
    if stypes['ies'] == 0:
        for i in range(stypes['n']):
            ti = stypes['vec_types'][i] - 1
            radius = stypes['radii'][ti] * stypes['frac_radii'][ti][stypes['num_layers'][ti] - 1]
            vtot += math.pow(radius, 3.0)
    else:
        radius = stypes['radii'][0] * stypes['frac_radii'][0][-1]
        vtot = math.pow(radius, 3.0)
    return (4.0 * math.pi / 3.0) * vtot

## \brief Parse the command line arguments.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'scat_name': '',
        'geom_name': '',
        'specific': [1.0],
        'medium': 0.0,
        'help_mode': False,
        'version_mode': False
    }
    skip_arg = False
    for arg_index in range(1, len(argv)):
        if skip_arg:
            skip_arg = False
            continue
        arg = argv[arg_index]
        split_arg = arg.split('=', 1)
        if (arg == "--scat" or arg == "--geom"):
            if (arg_index + 1 >= len(argv)):
                raise ValueError("%s needs a file name."%arg)
            config['scat_name' if arg == "--scat" else 'geom_name'] = argv[arg_index + 1]
            skip_arg = True
        elif (arg.startswith("--spew=")):
            try:
                config['specific'] = [float(value) for value in split_arg[1].split(',')]
            except ValueError:
                raise ValueError("invalid specific weights \'{0:s}\'".format(split_arg[1]))
        elif (arg.startswith("--medium=")):
            try:
                config['medium'] = float(split_arg[1])
            except ValueError:
                raise ValueError("invalid specific weight \'{0:s}\'".format(split_arg[1]))
        elif (arg == "--help"):
            config['help_mode'] = True
        elif (arg == "--version"):
            config['version_mode'] = True
        else:
            raise ValueError("Unrecognized argument \'{0:s}\'".format(arg))
    return config

## \brief Print a command-line help summary.
def print_help():
    print("                                                                           ")
    print("***               COMPUTE THE TENSOR OF INERTIA                          ***")
    print("                                                                           ")
    print("Compute the mass, the tensor of inertia and the filling factor of a        ")
    print("NP_TMcode particle model.                                                  ")
    print("                                                                           ")
    print("Usage: \"./inertia.py --scat SCAT --geom GEOM [OPTIONS]\"                  ")
    print("                                                                           ")
    print("Valid options are:                                                         ")
    print("--scat SCAT              Scatterer configuration file (mandatory).         ")
    print("--geom GEOM              Geometry configuration file (mandatory).          ")
    print("--spew=W1[,W2,...]       Specific weights of the materials, in g/cm3       ")
    print("                         (default 1). The materials are the distinct       ")
    print("                         dielectric constants of the first wavelength, in  ")
    print("                         the order of the scatterer file; one value applies")
    print("                         to all of them.                                   ")
    print("--medium=W               Specific weight of the external medium, which     ")
    print("                         fills the outer part of a last transition layer   ")
    print("                         (in g/cm3, default 0).                            ")
    print("--help                   Print this help and exit.                         ")
    print("--version                Print script version and exit.                    ")
    print("                                                                           ")

# ### PROGRAM EXECUTION ###
## \cond
if __name__ == "__main__":
    res = main()
    exit(res)
## \endcond
