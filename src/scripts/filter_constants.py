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

## @package filter_constants
#  \brief Script to create optical function grids for NPTM_code.
#
#  One of the factors that critically affect the duration of a calculation
#  with NPTM_code is the number of wavelengths for which the calculation
#  is executed. This number depends on the definition of the optical functions
#  of the involved materials. The filter_constants.py script can be used to
#  extract optimized grids of optical constants.
#
#  The script requires python3.

import math
import numpy as np
#import pdb
from sys import argv

## \cond
__version__ = "0.10.10"
allow_plots = True
## \endcond

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("WARNING: MATPLOTLIB not found - disabling plots.")
    allow_plots = False

## \brief Main execution code.
#
#  `main()` is the function that handles the creation of the code configuration.
#  It returns an integer value as exit code, using 0 to signal successful execution.
#
#  \returns result: `int` Exit code (0 = SUCCESS).
def main():
    result = 0
    config = parse_arguments()
    if config['version_mode']:
        print("filter_constants.py v%s."%__version__)
    elif (config['help_mode']):
        print_help()
    else:
        if (config['input_files'] == ''):
            print("ERROR: no input files provided (--in INPUT)!")
            result = 1
        else:
            split_input = config['input_files'].split(',')
            if (len(split_input) == 1):
                scan_info = scan_single_file(config, split_input[0])
                result = scan_info['exit_code']
                if (result == 0):
                    wl_factor = 1.0e6
                    wl_units = config['wl_units']
                    if (wl_units in ["nanometers", "nm"]):
                        wl_factor = 1.0e9
                    elif (wl_units in ["millimeters", "mm"]):
                        wl_factor = 1.0e3
                    elif (wl_units in ["centimeters", "cm"]):
                        wl_factor = 1.0e2
                    elif (wl_units in ["decimeters", "dm"]):
                        wl_factor = 1.0e1
                    elif (wl_units in ["meters", "m"]):
                        wl_factor = 1.0
                    write_output_file(config['output_files'], scan_info, wl_factor)
                    # INFO section
                    num_filtered_data = scan_info['num_filtered_data']
                    num_orig_data = scan_info['num_orig_data']
                    msg_filtered = "{0:d} filtered lines".format(num_filtered_data) if num_filtered_data != 1 else "1 filtered line"
                    msg_orig = "{0:d} input lines".format(num_orig_data) if num_orig_data != 1 else "1 input line"
                    max_distance = get_max_distance(scan_info)
                    msg_distance1 = "INFO: maximum absolute distance was {0:.5g} at {1:.5e} {2:s} in the {3:s} part".format(
                        max_distance['max_difference'],
                        max_distance['max_wavelength'],
                        config['wl_units'],
                        max_distance['differing_set']
                    )
                    msg_distance2 = "      (fitted value is {0:.5g}, actual data value is {1:.5g}).".format(
                        max_distance['max_fitted'],
                        max_distance['max_value']
                    )
                    print("INFO: extracted %s out of %s."%(msg_filtered, msg_orig))
                    print(msg_distance1)
                    print(msg_distance2)
                    if (config['make_plots']):
                        plot_data([scan_info], config['wl_units'])
                # end result == 0 check
            else:
                split_output = config['output_files'].split(',')
                aligned_infos = match_multiple_files(config)
                for i in range(len(aligned_infos)):
                    info = aligned_infos[i]
                    ecode = info['exit_code']
                    file_name = split_output[i]
                    if (ecode == 0):
                        write_output_file(file_name, info)
                        # INFO section
                        num_filtered_data = info['num_filtered_data']
                        num_orig_data = info['num_orig_data']
                        msg_filtered = "{0:d} filtered lines".format(num_filtered_data) if num_filtered_data != 1 else "1 filtered line"
                        msg_orig = "{0:d} input lines".format(num_orig_data) if num_orig_data != 1 else "1 input line"
                        max_distance = get_max_distance(info)
                        msg_distance1 = "INFO: maximum absolute distance was {0:.5g} at {1:.5e} {2:s} in the {3:s} part".format(
                            max_distance['max_difference'],
                            max_distance['max_wavelength'],
                            config['wl_units'],
                            max_distance['differing_set']
                        )
                        msg_distance2 = "      (fitted value is {0:.5g}, actual data value is {1:.5g}).".format(
                            max_distance['max_fitted'],
                            max_distance['max_value']
                        )
                        print("INFO: extracted %s out of %s."%(msg_filtered, msg_orig))
                        print(msg_distance1)
                        print(msg_distance2)
                        # end INFO section
                    else:
                        print("WARNING: scanning {0:s} resulted in error code {1:d}.".format(file_name, ecode))
                    result += ecode
                if (result == 0):
                    if (config['make_plots']):
                        plot_data(aligned_infos, config['wl_units'])
    return result

## \brief Maximum distance between input data and linear interpolation of filtered data.
#  
#  This function returns a dictionary object to make a quick estimate of the quality
#  of the chosen filter. The diagnostic is the maximum absolute offset between the
#  distribution of the input data and a segmented linear interpolation touching all
#  the filtered data.
#
#  \param[in] scan_info: `dict` A dictionary containing the results of a file scan.
#  \returns result: `dict` A dictionary containing the wavelength of the maximum difference
#           (`max_wavelength`), the value of the input data at that wavelength (`mav_value`),
#           the value of the interpolated filtered functions at the same wavelength (`max_fitted`),
#           the offset between interpolation and data (`max_difference`), and the set of values
#           where the difference was observed (`differing_set`, being either REAL or IMAGINARY).
def get_max_distance(scan_info):
    wl_orig = scan_info['wl_orig']
    reps_orig = scan_info['reps_orig']
    ieps_orig = scan_info['ieps_orig']
    wl_filtered = scan_info['wl_filtered']
    reps_filtered = scan_info['reps_filtered']
    ieps_filtered = scan_info['ieps_filtered']
    max_difference = 0.0
    max_fitted = 0.0
    max_value = 0.0
    max_wavelength = 0.0
    differing_set = ""
    filtered_index = 1
    for i in range(len(wl_orig)):
        wl = wl_orig[i]
        if (wl <= wl_filtered[0]):
            continue
        if (wl < wl_filtered[filtered_index]):
            dx = wl - wl_filtered[filtered_index - 1]
            dry = reps_filtered[filtered_index] - reps_filtered[filtered_index - 1]
            diy = ieps_filtered[filtered_index] - ieps_filtered[filtered_index - 1]
            dwl = wl_filtered[filtered_index] - wl_filtered[filtered_index - 1]
            rp = reps_filtered[filtered_index - 1] + dry * dx / dwl
            ip = ieps_filtered[filtered_index - 1] + diy * dx / dwl
            rdiff = rp - reps_orig[i]
            idiff = ip - ieps_orig[i]
            if (rdiff < 0.0):
                rdiff *= -1.0
            if (idiff < 0.0):
                idiff *= -1.0
            if (rdiff > max_difference):
                max_difference = rdiff
                max_fitted = rp
                max_value = reps_orig[i]
                max_wavelength = wl
                differing_set = "real"
            if (idiff > max_difference):
                max_difference = idiff
                max_fitted = ip
                max_value = ieps_orig[i]
                max_wavelength = wl
                differing_set = "imaginary"
        else:
            filtered_index += 1
            if (filtered_index == len(wl_filtered)):
                break # for i
        # end of wl step check
    result = {
        'max_wavelength': max_wavelength,
        'max_difference': max_difference,
        'max_fitted': max_fitted,
        'max_value': max_value,
        'differing_set': differing_set
    }
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
        'input_files': '',
        'output_files': '',
        'help_mode': False,
        'version_mode': False,
        'force_peaks': True,
        'make_plots': allow_plots,
        'step': 5.0e-8,
        'threshold': 0.1,
        'wl_start': 0.0,
        'wl_end': 0.0,
        'wl_tolerance': 0.0,
        'wl_units': 'micrometers'
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
        elif (arg.startswith("--no-peaks")):
            config['force_peaks'] = False
        elif (arg.startswith("--no-plots")):
            config['make_plots'] = False
        elif (arg.startswith("--in")):
            dict_key = 'input_files'
            skip_arg = True
        elif (arg.startswith("--out")):
            dict_key = 'output_files'
            skip_arg = True
        elif (arg.startswith("--step=")):
            split_arg = arg.split('=')
            config['step'] = float(split_arg[1])
        elif (arg.startswith("--threshold=")):
            split_arg = arg.split('=')
            config['threshold'] = float(split_arg[1])
        elif (arg.startswith("--wl-start=")):
            split_arg = arg.split('=')
            config['wl_start'] = float(split_arg[1])
        elif (arg.startswith("--wl-stop=")):
            split_arg = arg.split('=')
            config['wl_end'] = float(split_arg[1])
        elif (arg.startswith("--wl-tol=")):
            split_arg = arg.split('=')
            config['wl_tolerance'] = float(split_arg[1])
        elif (arg.startswith("--wl-units=")):
            known_units = [
                'nanometers', 'nm', 'micrometers', 'um', 'millimeters', 'mm',
                'centimeters', 'cm', 'decimeters', 'dm', 'meters', 'm'
            ]
            split_arg = arg.split('=')
            config['wl_units'] = split_arg[1]
            if (config['wl_units'] not in known_units):
                raise Exception("Unrecognized wavelength units %s!"%config['wl_units'])
        else:
            raise Exception("Unrecognized argument \"{0:s}\"!".format(arg))
    # end for loop
    if (config['output_files'] == ''):
        split_in = config['input_files'].split(',')
        for name_in in split_in:
            name_out = name_in.split('.')[0] + "_filtered.csv"
            if (config['output_files'] == ''):
                config['output_files'] = name_out
            else:
                config['output_files'] += ",{0:s}".format(name_out)
    else:
        split_in = config['input_files'].split(',')
        split_out = config['output_files'].split(',')
        if (len(split_in) != len(split_out)):
            raise Exception("Output list does not match input!")
        for ni in range(len(split_in)):
            input_name = split_in[ni]
            output_name = split_out[ni]
            if (input_name == output_name):
                raise Exception("No input file overwriting allowed!")
    return config

## \brief Make a quick-look plot with MATPLOTLIB.
#
#  \param[in] scan_infos: `list` A list of dictionaries with the scanned data file.
#  \param[in] wl_units: `string` The name of the units for the wavelength scale.
def plot_data(scan_infos, wl_units):
    cname_array = [
        'red',
        'blue',
        'green',
        'purple',
        'orange',
        'cyan',
        'grey',
        'black'
    ]
    # Plot making section
    for i in range(len(scan_infos)):
        scan_info = scan_infos[i]
        wl_orig = scan_info['wl_orig']
        reps_orig = scan_info['reps_orig']
        ieps_orig = scan_info['ieps_orig']
        wl_filtered = scan_info['wl_filtered']
        reps_filtered = scan_info['reps_filtered']
        ieps_filtered = scan_info['ieps_filtered']
        rcname = cname_array[(2 * i) % 8]
        icname = cname_array[(2 * i + 1) % 8]
        plt.plot(wl_orig, reps_orig, color=rcname, marker='', ls='-', label=r"Original $\mathfrak{Re}(\varepsilon)$")
        plt.plot(wl_filtered, reps_filtered, color=rcname, marker='o', ls='', label=r"Filtered $\mathfrak{Re}(\varepsilon)$")
        plt.plot(wl_orig, ieps_orig, color=icname, marker='', ls='--', label=r"Original $\mathfrak{Im}(\varepsilon)$")
        plt.plot(wl_filtered, ieps_filtered, color=icname, marker='s', ls='', label=r"Filtered $\mathfrak{Im}(\varepsilon)$")
    plt.xlabel("Wavelength ({0:s})".format(wl_units))
    plt.ylabel(r"$\mathfrak{Re}(\varepsilon)$|$\mathfrak{Im}(\varepsilon)$")
    plt.legend(loc="best")
    plt.show()
    # end plot making section

## \brief Print a command-line help summary.
def print_help():
    print("      ###############################################                 ")
    print("      #                                             #                 ")
    print("      #          NPTM_code FILTER_CONSTANTS         #                 ")
    print("      #                                             #                 ")
    print("      ###############################################                 ")
    print("                                                                      ")
    print("Filter the optical functions for a model material.                    ")
    print("                                                                      ")
    print("Usage: \"./filter_constants.py --in INPUT [options]\"                 ")
    print("                                                                      ")
    print("Valid options are:                                                    ")
    print("--in                  Comma separated list of input files (mandatory).")
    print("--out                 Comma separated list of output files (optional, ")
    print("                      but, if given, must be as many as for --in).    ")
    print("--help                Print this help and exit.                       ")
    print("--no-peaks            Disable mandatory extraction of peaks (default  ")
    print("                      is enabling peaks).                             ")
    print("--no-plots            Disable plotting the solution with MATPLOTLIB   ")
    print("                      (default is enable plots).                      ")
    print("--step=VALUE          Regular step in meters to sample flat regions   ")
    print("                      of the functions (use <= 0 to disable; default  ")
    print("                      is 5e-8).                                       ")
    print("--threshold=VALUE     Relative tolerance threshold to pick a point in ")
    print("                      the filtering process (optional, default is     ")
    print("                      0.1).                                           ")
    print("--wl_start=VALUE      Starting wavelength in meters for the filtering ")
    print("                      window.                                         ")
    print("--wl_stop=VALUE       Ending wavelength in meters for the filtering   ")
    print("                      window.                                         ")
    print("--wl_tol=VALUE        Minimum separation to consider two wavelength   ")
    print("                      values as distinct in multiple files (default is")
    print("                      0.001 times the shortest wavelength).")
    print("--wl_units=UNITS      Name of the wavelength units ONLY FOR PLOTTING  ")
    print("                      PURPOSES (data must be always in meters, only   ")
    print("                      MATPLOTLIB uses this setting for formatting).   ")
    print("--version             Print script version and exit.                  ")
    print("                                                                      ")

## \brief Filter a single file based on the configuration options.
#
#  Perform filtering of a single file based on custom thresholda and step
#  configurations. The filtered data are written to a CSV file, then a
#  diagnostic log is printed to terminal. Optionally the filter selection
#  is shown as a plot, if MATPLOTLIB is available on the system.
#
#  \param[in] config: `dict` A dictionary containing the script configuration.
#  \param[in] file_name: `string` The name of the single input file.
#  \return result: `dict` A dictionary containing the results of the scan,
#          including "exit_code" (`int`, 0 if succesful), "wl_orig" (`array-like`,
#          the original wavelength scale), "reps_orig" (`array-like`, the original
#          real parts of the dielectric functions), "ieps_orig" (`array-like`, the
#          original imaginary parts of the dielectric functions), "wl_filtered"
#          (`array-like`, the filtered wavelength scale), "reps_filtered"
#          (`array-like`, the filtered real parts of the dielectric functions),
#          "ieps_filtered" (`array-like`, the filtered imaginary parts of the
#          dielectric functions), and "reason" (`array-like`, containing a code
#          to track whether a point was collected for step reasons [1], for
#          threshold filter [2], or for being a peak point [3]).
def scan_single_file(config, file_name):
    result = {
        'exit_code': -1,
        'header': "",
        'num_read_lines': 0,
        'num_orig_data': 0,
        'num_filtered_data': 0,
        'wl_orig': [],
        'reps_orig': [],
        'ieps_orig': [],
        'wl_filtered': [],
        'reps_filtered': [],
        'ieps_filtered': [],
        'reason': []
    }
    try:
        input_file = open(file_name, 'r')
        file_line = input_file.readline()
        num_read_lines = 1
        num_orig_data = 0
        num_filtered_data = 0
        wl_factor = 1.0e6
        wl_units = config['wl_units']
        if (wl_units in ["nanometers", "nm"]):
            wl_factor = 1.0e9
        elif (wl_units in ["millimeters", "mm"]):
            wl_factor = 1.0e3
        elif (wl_units in ["centimeters", "cm"]):
            wl_factor = 1.0e2
        elif (wl_units in ["decimeters", "dm"]):
            wl_factor = 1.0e1
        elif (wl_units in ["meters", "m"]):
            wl_factor = 1.0
        wl_orig = result['wl_orig']
        reps_orig = result['reps_orig']
        ieps_orig = result['ieps_orig']
        wl_filtered = result['wl_filtered']
        reps_filtered = result['reps_filtered']
        ieps_filtered = result['ieps_filtered']
        reason = result['reason']
        step = config['step']
        threshold = config['threshold']
        wl0 = 0.0
        wl1 = 0.0
        reps0 = 0.0
        reps1 = 0.0
        ieps0 = 0.0
        ieps1 = 0.0
        last_dreps = 0.0
        last_dieps = 0.0
        can_write = True
        while (file_line != ""):
            if (file_line.startswith('#')):
                result['header'] += file_line
                file_line = input_file.readline()
                num_read_lines += 1
                continue
            split_line = file_line.split(',')
            if (len(split_line) == 3):
                # parse the line
                if (wl0 == 0.0):
                    # always parse the first data line
                    wl0 = float(split_line[0])
                    reps0 = float(split_line[1])
                    ieps0 = float(split_line[2])
                    wl_orig.append(wl0 * wl_factor)
                    reps_orig.append(reps0)
                    ieps_orig.append(ieps0)
                    num_orig_data += 1
                    if (wl0 >= config['wl_start']):
                        wl_filtered.append(wl0 * wl_factor)
                        reps_filtered.append(reps0)
                        ieps_filtered.append(ieps0)
                        reason.append(3)
                        num_filtered_data += 1
                else:
                    can_write = True
                    wl1 = float(split_line[0])
                    reps1 = float(split_line[1])
                    ieps1 = float(split_line[2])
                    dreps = reps1 - reps0
                    dieps = ieps1 - ieps0
                    wl_orig.append(wl1 * wl_factor)
                    reps_orig.append(reps1)
                    ieps_orig.append(ieps1)
                    num_orig_data += 1
                    if (wl1 < config['wl_start']):
                        num_read_lines += 1
                        file_line = input_file.readline()
                        continue # while loop
                    if (config['wl_end'] > config['wl_start'] and wl1 > config['wl_end']):
                        break # while loop
                    if (step > 0.0):
                        if (wl1 - wl0 >= step):
                            # compute the values at step location with linear interpolation
                            wl = wl0 + step
                            x0 = wl_orig[-2] / wl_factor if len(wl_orig) > 1 else wl0
                            x1 = wl_orig[-1] / wl_factor if len(wl_orig) > 1 else wl1
                            dx = wl - x0
                            ry0 = reps_orig[-2] if len(reps_orig) > 1 else reps0
                            ry1 = reps_orig[-1] if len(reps_orig) > 1 else reps1
                            dry = ry1 - ry0
                            iy0 = ieps_orig[-2] if len(ieps_orig) > 1 else ieps0
                            iy1 = ieps_orig[-1] if len(ieps_orig) > 1 else ieps1
                            diy = iy1 - iy0
                            reps = ry0 + dry * dx / (x1 - x0)
                            ieps = iy0 + diy * dx / (x1 - x0)
                            # write a line if step is enabled and satisfied
                            wl_filtered.append(wl * wl_factor)
                            reps_filtered.append(reps)
                            ieps_filtered.append(ieps)
                            reason.append(1)
                            num_filtered_data += 1
                            can_write = False
                            wl0 = wl
                            reps0 = reps
                            ieps0 = ieps
                    # end of step > 0.0 check
                    if (config['force_peaks']):
                        rpeak = (dreps * last_dreps < 0.0)
                        ipeak = (dieps * last_dieps < 0.0)
                        if ((rpeak or ipeak) and can_write):
                            # write a line if peaks are enabled and satisfied
                            can_write = False
                            wl_filtered.append(wl1 * wl_factor)
                            reps_filtered.append(reps1)
                            ieps_filtered.append(ieps1)
                            reason.append(3)
                            num_filtered_data += 1
                            wl0 = wl1
                            reps0 = reps1
                            ieps0 = ieps1
                    # end of force_peaks check
                    last_dreps = dreps
                    last_dieps = dieps
                    if (reps0 != 0.0):
                        rel_dreps = (reps0 + dreps) / reps0 if dreps > 0.0 else (reps0 - dreps) / reps0
                        if (rel_dreps < 0.0):
                            rel_dreps *= -1.0
                        if ((rel_dreps > 1.0 + threshold or rel_dreps < 1.0 - threshold) and can_write):
                            # write a line if tolerance is violated
                            can_write = False
                            wl_filtered.append(wl1 * wl_factor)
                            reps_filtered.append(reps1)
                            ieps_filtered.append(ieps1)
                            reason.append(2)
                            num_filtered_data += 1
                            wl0 = wl1
                            reps0 = reps1
                            ieps0 = ieps1
                    # end of reps0 != 0.0 check
                    if (ieps0 != 0.0):
                        rel_dieps = (ieps0 + dieps) / ieps0 if dieps > 0.0 else (ieps0 - dieps) / ieps0
                        if (rel_dieps < 0.0):
                            rel_dieps *= -1.0
                        if ((rel_dieps > 1.0 + threshold or rel_dieps < 1.0 - threshold) and can_write):
                            # write a line if tolerance is violated
                            can_write = False
                            wl_filtered.append(wl1 * wl_factor)
                            reps_filtered.append(reps1)
                            ieps_filtered.append(ieps1)
                            reason.append(2)
                            num_filtered_data += 1
                            wl0 = wl1
                            reps0 = reps1
                            ieps0 = ieps1
                    # end of reps0 != 0.0 check
                # end of wl0 == 0.0 check
            else:
                print("ERROR: invalid input file %s at line %d!"%(file_name, num_read_lines))
                result['exit_code'] = 1
                break # while loop
            # end of len(split_line) check
            file_line = input_file.readline()
            num_read_lines += 1
        # end of while loop
        if (wl1 <= config['wl_end'] and can_write):
            can_write = False
            wl_filtered.append(wl1 * wl_factor)
            reps_filtered.append(reps1)
            ieps_filtered.append(ieps1)
            reason.append(3)
            num_filtered_data += 1
        input_file.close()
        if (result['exit_code'] < 0):
            reason[-1] = 3
            result['exit_code'] = 0
            result['num_read_lines'] = num_read_lines
            result['num_orig_data'] = num_orig_data
            result['num_filtered_data'] = num_filtered_data
    except FileNotFoundError as ex:
        print("ERROR: file not found %s!"%config['input_files'])
        result['exit_code'] = 1
    return result

## \brief Filter multiple files based on the configuration options.
#
#  A sequence of optical function data files is filtered according to the
#  same samplig grid, resulting in a set of files ready for use in the same
#  simulation.
#
#  \param[in] config: `dict` A dictionary containing the script configuration.
#  \return aligned_infos: `list` A list of dictionaries containing the results
#          of filtering aligned to a common scale.
def match_multiple_files(config):
    scan_infos = []
    split_input = config['input_files'].split(',')
    wl_factor = 1.0e6
    wl_units = config['wl_units']
    if (wl_units in ["nanometers", "nm"]):
        wl_factor = 1.0e9
    elif (wl_units in ["millimeters", "mm"]):
        wl_factor = 1.0e3
    elif (wl_units in ["centimeters", "cm"]):
        wl_factor = 1.0e2
    elif (wl_units in ["decimeters", "dm"]):
        wl_factor = 1.0e1
    elif (wl_units in ["meters", "m"]):
        wl_factor = 1.0
    for file_name in split_input:
        scan_infos.append(scan_single_file(config, file_name))
    # end scan_infos loop
    # Find the global X range
    x_min = min(np.min(s['wl_filtered']) for s in scan_infos)
    x_max = max(np.max(s['wl_filtered']) for s in scan_infos)
    tolerance = config['wl_tolerance'] if config['wl_tolerance'] != 0.0 else (1.0e-3 * x_min / wl_factor)
    # Extract special points
    special_ieps = []
    for s in scan_infos:
        x_arr = np.array(s['wl_filtered'])
        info_arr = np.array(s['reason'])
        special_ieps.extend(x_arr[info_arr > 1])
    # Make a regular grid
    num_regular_points = int((x_max - x_min) / (config['step'] * wl_factor))
    regular_x = np.linspace(x_min, x_max, num_regular_points)
    # Get a coarse global X vector
    coarse_x = np.sort(np.unique(np.concatenate([regular_x, special_ieps])))
    # Numerical tolerance filtering
    mask = np.insert(np.diff(coarse_x) > tolerance, 0, True)
    common_x = coarse_x[mask]
    # Re-align each series on the common scale
    aligned_infos = []
    for s in scan_infos:
        wl_old = np.array(s['wl_filtered'])
        reps_old = np.array(s['reps_filtered'])
        ieps_old = np.array(s['ieps_filtered'])
        info_old = np.array(s['reason'])
        
        # Value interpolation on the new grid
        reps_interp = np.interp(common_x, wl_old, reps_old)
        ieps_interp = np.interp(common_x, wl_old, ieps_old)
        
        # Mapping of original INFO on the new scale
        info_new = np.ones(len(common_x), dtype=int)
        for x_val, info_val in zip(wl_old, info_old):
            if info_val > 1:
                # Find corresponding index in common_x
                idx = np.argmin(np.abs(common_x - x_val))
                if np.abs(common_x[idx] - x_val) <= tolerance:
                    info_new[idx] = info_val
        # end of x_val, info_val loop

        aligned_infos.append({
            'exit_code': 0,
            'header': s['header'],
            'num_read_lines': s['num_read_lines'],
            'num_orig_data': s['num_orig_data'],
            'num_filtered_data': s['num_filtered_data'],
            'wl_orig': np.array(s['wl_orig']),
            'reps_orig': np.array(s['reps_orig']),
            'ieps_orig': np.array(s['ieps_orig']),
            'wl_filtered': common_x,
            'reps_filtered': reps_interp,
            'ieps_filtered': ieps_interp,
            'reason': info_new
        })
    # end of scan_infos loop
    return aligned_infos

## \brief Write the filtered data to an output file.
#
#  \param[in] file_name: `string` The name of the file to be written.
#  \param[in] scan_info: `dict` A dictionary with the scanned data file.
#  \param[in] wl_factor: `float` Wavelength unit conversion to meters.
def write_output_file(file_name, scan_info, wl_factor):
    output_file = open(file_name, 'w')
    output_file.write(scan_info['header'])
    for i in range(len(scan_info['wl_filtered'])):
        file_line = "{0:.5E},{1:.5E},{2:.5E}\n".format(
            scan_info['wl_filtered'][i] / wl_factor,
            scan_info['reps_filtered'][i],
            scan_info['ieps_filtered'][i]
        )
        output_file.write(file_line)
    output_file.close()
    
## \brief Exit code (0 for success).
exit_code = main()
exit(exit_code)
