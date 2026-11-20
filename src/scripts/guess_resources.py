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
#  \ingroup python_tools
#  \brief Script to estimate the computational cost of a model.
#
#  Script to assist in the creation of model input files starting
#  from a YAML descriptor.
#
#  The script requires python3.

import math
import numpy as np
import re
from sys import argv

import nptm_estimates as est

## \cond
from nptm_version import __version__
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
                print_resource_summary(model, config)
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
        'version_mode': False,
        'backend': "lapack",
        'threads': 1,
        'processes': 1,
        'gpus': 1,
        'host_ram': 0.0,
        'gpu_ram': 0.0,
        'no_lapack': False,
        'int_bytes': 8
    }
    value_keys = {
        "--edfb": 'edfb_file', "--geom": 'geom_file', "--backend": 'backend',
        "--threads": 'threads', "--processes": 'processes', "--gpus": 'gpus',
        "--host-ram": 'host_ram', "--gpu-ram": 'gpu_ram'
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
        elif (arg == "--no-lapack"):
            config['no_lapack'] = True
        elif (arg == "--lp64"):
            config['int_bytes'] = 4
        elif (arg in value_keys):
            dict_key = value_keys[arg]
            skip_arg = True
        else:
            raise Exception("Unrecognized argument \'{0:s}\'".format(arg))
    for key in ('threads', 'processes', 'gpus'):
        config[key] = int(config[key])
    for key in ('host_ram', 'gpu_ram'):
        config[key] = float(config[key])
    if (config['backend'] not in est.BACKENDS):
        raise Exception("Unknown back end \'{0:s}\' (use one of {1:s})".format(
            config['backend'], ", ".join(est.BACKENDS)))
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
    print("Usage: \"./guess_resources.py --edfb EDFB --geom GEOM [OPTIONS]\"")
    print("                                                          ")
    print("EDFB and GEOM must be valid NPTM_code configuration files.")
    print("The memory estimates are those that the programs log and  ")
    print("check; the runtime options of GEOM (USE_REFINEMENT,       ")
    print("INVERSION, DEVICE_POLICY, HOST_RAM_GB, GPU_RAM_GB) apply. ")
    print("                                                          ")
    print("Valid options are:                                        ")
    print("--backend BACKEND     Inversion back end of the build:    ")
    print("                      magma, cublas, lapack (default) or  ")
    print("                      lucin (no LAPACK).                  ")
    print("--threads N           OpenMP threads per MPI process      ")
    print("                      (default 1).                        ")
    print("--processes N         MPI processes per node (default 1). ")
    print("--gpus N              GPUs per node (default 1).          ")
    print("--host-ram GIB        Host memory of a node, to check the ")
    print("                      estimate (HOST_RAM_GB overrides it).")
    print("--gpu-ram GIB         Free memory of a GPU, to check the  ")
    print("                      estimate (GPU_RAM_GB overrides it). ")
    print("--no-lapack           A MAGMA or cuBLAS build without     ")
    print("                      LAPACK (CPU fall-back with lucin).  ")
    print("--lp64                32-bit LAPACK and MAGMA integers.   ")
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
    # As in the code: the smallest sphere centred in the origin of the
    # coordinates that contains the particle (the external sphere of an
    # inclusion), which sets the external order with dynamic orders.
    radii = [model['ros'][model['vec_types'][i] - 1] for i in range(model['nsph'])]
    Rcirc = est.particle_radius(
        radii, model['vec_x'], model['vec_y'], model['vec_z'], model['ext_radius']
    )
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
    print("INFO: minimum encircling radius around the origin is Rcirc = %.5em"%Rcirc)
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

## \brief Format a size in bytes as GiB, as the code logs it.
#
#  \param size: `int` Size in bytes.
#  \return text: `str` The size in GiB, with 5 significant digits.
def gib(size):
    return "{0:.5g} GiB".format(size / est.GIB)

## \brief Print a summary of the memory needed by the calculation.
#
#  The estimates are those that the programs log and check before the
#  calculation (see `nptm_estimates`). The host memory of a node is the sum
#  over its MPI processes; process 0 also holds the results of all the
#  scales. Each GPU holds the inversion of the threads of the processes that
#  use it, plus `GPU_PROCESS_OVERHEAD` per process, and the programs allow
#  `GPU_USABLE_FRACTION` of its free memory.
#
#  \param[in] model: `dict` A model description dictionary.
#  \param[in] config: `dict` The script configuration (build and run settings).
def print_resource_summary(model, config):
    threads = config['threads']
    processes = config['processes']
    gpus = config['gpus']
    backend = config['backend']
    needs = est.memory_needs(
        model, backend, threads, not config['no_lapack'], config['int_bytes']
    )
    print("INFO: estimates for a %s build, %d MPI process(es) per node, %d thread(s) each."%(
        backend, processes, threads))
    if (model['app'] != "SPHERE"):
        print("INFO: refinement is %s; inversion mode is %s."%(
            "on" if model['refine'] else "off", model['inv_mode']))
        print("INFO: size of matrices to invert: %d x %d (%s)."%(
            needs['matrix_order'], needs['matrix_order'], gib(needs['matrix'])))
    print("INFO: iteration data require %s of RAM per thread."%gib(needs['iteration_data']))
    if (model['app'] != "SPHERE"):
        print("INFO: matrix inversion requires %s of RAM per thread."%gib(needs['workspace']))
    print("INFO: output data require %s of RAM."%gib(needs['output']))
    print("INFO: MPI process 0 needs %s of RAM (%d thread(s))."%(gib(needs['process_0']), threads))
    node_need = needs['process_0'] + (processes - 1) * needs['other_process']
    if (processes > 1):
        print("INFO: each other MPI process needs %s of RAM."%gib(needs['other_process']))
        print("INFO: the node of MPI process 0 needs %s of RAM."%gib(node_need))
    host_ram = model['host_ram_gb'] if model['host_ram_gb'] > 0.0 else config['host_ram']
    if (host_ram > 0.0):
        if (node_need <= host_ram * est.GIB):
            print("INFO: it fits in %.5g GiB of host memory."%host_ram)
        else:
            print("WARNING: it does not fit in %.5g GiB of host memory."%host_ram)
    if (needs['device_per_thread'] > 0 and gpus > 0):
        sharing = (processes + gpus - 1) // gpus
        gpu_need = sharing * needs['device_per_process']
        print("INFO: each thread needs %s on the GPU for the inversion."%gib(needs['device_per_thread']))
        print("INFO: the calculation needs %s of memory on each GPU (with %s for the GPU runtime)."%(
            gib(gpu_need), gib(sharing * est.GPU_PROCESS_OVERHEAD)))
        if (model['gpu_ram_gb'] > 0.0):
            usable = model['gpu_ram_gb'] * est.GIB
            source = "allowed by GPU_RAM_GB"
        elif (config['gpu_ram'] > 0.0):
            usable = est.GPU_USABLE_FRACTION * config['gpu_ram'] * est.GIB
            source = "usable ({0:d}% of {1:.5g} GiB free)".format(
                int(100 * est.GPU_USABLE_FRACTION + 0.5), config['gpu_ram'])
        else:
            usable = 0.0
            print("INFO: it needs at least %s of free GPU memory (%d%% usable)."%(
                gib(gpu_need / est.GPU_USABLE_FRACTION), int(100 * est.GPU_USABLE_FRACTION + 0.5)))
        if (usable > 0.0):
            if (gpu_need <= usable):
                print("INFO: it fits in %s %s."%(gib(usable), source))
            else:
                print("WARNING: it does not fit in %s %s."%(gib(usable), source))

## \brief Scan the calculation model.
#
#  The computational costs depennd the characteristics of the model. This
#  function scans the model configuration files and returns a model description
#  dictionary.
#
#  \param edfb_name: `str` Path to the scatterer configuration (EDFB) file.
#  \param geom_name: `str` Path to the geometry configuration file.
#  \returns model: `dict` A dictionary containing the model description.
def scan_model(edfb_name, geom_name):
    file_line = "INIT"
    model = {
        'app': "",
        'nsph': 0,
        'ies': 0,
        'lm': 0,
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
        'vec_z': [],
        'ext_radius': None,
        'angles': [],
        'isam': 0,
        'refine': False,
        'inv_mode': "LU",
        'device_policy': "STOP",
        'host_ram_gb': 0.0,
        'gpu_ram_gb': 0.0
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
            # one line with the first value and the step
            edfb_file.readline()
            read_lines += 1
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
                # the code expands the external field on the external sphere
                model['ext_radius'] = radius
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
                # one line with the first value and the step
                edfb_file.readline()
                read_lines += 1
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
                # one line with the first value and the step
                edfb_file.readline()
                read_lines += 1
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
        model['isam'] = int(array_values[8].group())
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
        model['lm'] = model['li']
        model['npnt'] = int(array_values[3].group())
        model['npntts'] = int(array_values[4].group())
        model['isam'] = int(array_values[5].group())
        model['vec_x'].append(0.0)
        model['vec_y'].append(0.0)
        model['vec_z'].append(0.0)
    # Incident and scattered theta, then incident and scattered phi.
    for ai in range(2):
        file_line = geom_file.readline()
        read_lines += 1
        str_line = file_line.replace('D', 'E').replace('d', 'E')
        array_values = iter_to_array(number_reg.finditer(str_line))
        if (len(array_values) != 6):
            print("ERROR: %s is not a valid geometry configuration file!"%geom_name)
            print("   INVALID LINE: \"%s\""%file_line[:-1])
            print("   at line %d in %s."%(read_lines, geom_name))
            geom_file.close()
            return None
        for av in array_values:
            model['angles'].append(float(av.group()))
    # Runtime options (KEY=VALUE lines after the mandatory data).
    for file_line in geom_file:
        option = file_line.strip()
        if (option.startswith("USE_REFINEMENT=")):
            model['refine'] = (int(option.split("=", 1)[1]) > 0)
        elif (option.startswith("INVERSION=")):
            model['inv_mode'] = option.split("=", 1)[1].strip()
        elif (option.startswith("DEVICE_POLICY=")):
            model['device_policy'] = option.split("=", 1)[1].strip()
        elif (option.startswith("HOST_RAM_GB=")):
            model['host_ram_gb'] = float(option.split("=", 1)[1])
        elif (option.startswith("GPU_RAM_GB=")):
            model['gpu_ram_gb'] = float(option.split("=", 1)[1])
    geom_file.close()
    # Names used by nptm_estimates.
    model['nshl'] = model['vec_nshl']
    model['configurations'] = len(model['vec_nshl'])
    return model

## \brief Exit code (0 for success).
exit_code = main()
exit(exit_code)
