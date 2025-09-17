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

## @package pydynrange
#  \brief Script to calculate the dynamic range of a complex matrix.
#
#  The script execution requires python3.

import cmath
import math
import os

from sys import argv, stdout

## \brief Class to represent the dynamic range memory structure.
#
# In order to implement a light-weight test that is able to choose
# whether to produce a text-only output or a text + diagnostic plot
# output, the `PlotData` class stores the memory structures to
# produce the text log and the plots. The data for the plots are
# filled with dummy place-holders if the required output is text
# only, while they are filled with actual data, otherwise.
#
# \returns exit_code: `int` 0 on successful completion.
class PlotData:
    ## \brief PlotData instance constructor.
    def __init__(self):
        ## \brief Maximum real part value.
        self.max_real = -1.0
        ## \brief Maximum imaginary part value.
        self.max_imag = -1.0
        ## \brief Maximum absolute value.
        self.max_absl = -1.0
        ## \brief Minimum (most negative) real part value.
        self.min_real = 1.0
        ## \brief Minimum (most negative) imaginary part value.
        self.min_imag = 1.0
        ## \brief Smallest absolute value.
        self.sml_absl = 1.0
        ## \brief Smallest absolute real part value.
        self.sml_real = 1.0
        ## \brief Smallest absolute imaginary part value.
        self.sml_imag = 1.0
        ## \brief Histogram of absolute values.
        self.absl_hist = [0 for i in range(-99, 100)]
        ## \brief Histogram of real values.
        self.real_hist = [0 for i in range(-99, 100)]
        ## \brief Histogram of imaginary values.
        self.imag_hist = [0 for i in range(-99, 100)]
        ## \brief Maximum real order of magnitude.
        self.max_real_mag = -100
        ## \brief Minimum real order of magnitude.
        self.min_real_mag = 100
        ## \brief Maximum imaginary order of magnitude.
        self.max_imag_mag = -100
        ## \brief Minimum imaginary order of magnitude.
        self.min_imag_mag = 100
        ## \brief Maximum absolute order of magnitude.
        self.max_absl_mag = -100
        ## \brief Minimum absolute order of magnitude.
        self.min_absl_mag = 100

    ## \brief Print a text log of the dynamic range.
    def log_dynamic_range(self):
        print("        MAX( ABS[AM] ) = %14.7e"%self.max_absl)
        print("        MIN( ABS[AM] ) = %14.7e"%self.sml_absl)
        print("       MAX( REAL[AM] ) = %14.7e"%self.max_real)
        print("       MIN( REAL[AM] ) = %14.7e"%self.min_real)
        print("       MAX( IMAG[AM] ) = %14.7e"%self.max_imag)
        print("       MIN( IMAG[AM] ) = %14.7e"%self.min_imag)
        print("MIN( ABS( REAL[AM] ) ) = %14.7e"%self.sml_real)
        print("MIN( ABS( IMAG[AM] ) ) = %14.7e"%self.sml_imag)
        return

    ## \brief Make histograms of dynamic range with matplotlib.
    #
    #  \param config: `dict` Dictionary of configuration options.
    def plot_dynamic_range(self, config):
        import matplotlib.pyplot as plt
        import numpy as np
        self.absl_hist = np.array(self.absl_hist)
        self.real_hist = np.array(self.real_hist)
        self.imag_hist = np.array(self.imag_hist)
        vec_x_pos = np.arange(-98.5, 100.0, 1.0)
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)
        ax1.bar(vec_x_pos, self.absl_hist, color = "grey", label = "Abs. values")
        ax1.set_xlim(self.min_absl_mag - 1, self.max_absl_mag + 1)
        ax1.set_xlabel("Order of magnitude")
        ax1.set_ylabel("Element counts")
        ax1.legend(loc = "best")
        ax2.bar(vec_x_pos, self.real_hist, color = "red", label = "Real values")
        ax2.set_xlim(self.min_real_mag - 1, self.max_real_mag + 1)
        ax2.set_xlabel("Order of magnitude")
        ax2.set_ylabel("Element counts")
        ax2.legend(loc = "best")
        ax3.bar(vec_x_pos, self.imag_hist, color = "blue", label = "Imag. values")
        ax3.set_xlim(self.min_imag_mag - 1, self.max_imag_mag + 1)
        ax3.set_xlabel("Order of magnitude")
        ax3.set_ylabel("Element counts")
        ax3.legend(loc = "best")
        if (config['log_plot_y']):
            ax1.set_yscale("log")
            ax2.set_yscale("log")
            ax3.set_yscale("log")
        plt.tight_layout()
        plt.show()
        return

## \brief Main execution code
#
# `main()` is the function that handles the creation of the script configuration
# and the execution of the calculation. It returns 0 on successful completion.
#
# \returns exit_code: `int` 0 on successful completion.
def main():
    config = parse_arguments()
    exit_code = 0
    if config['help_mode'] or len(argv) == 1:
        config['help_mode'] = True
        print_help()
    else:
        if config['matrix_name'] == "":
            exit_code = 1
        else:
            exit_code = get_dynamic_range(config)
    return exit_code

## \brief Compute the dynamic range of a matrix.
#
#  \param config: `dict` A dictionary containing the script configuration.
#
#  \returns result: `int` The exit code of the operation (0 for success).
def get_dynamic_range(config):
    result = 0
    num_read_lines = 0
    pd = PlotData()
    stdout.write("INFO: scanning data file. Please, wait...")
    stdout.flush()
    matrix_file = open(config['matrix_name'], 'r')
    str_line = matrix_file.readline()
    str_line = matrix_file.readline()
    str_line = matrix_file.readline()
    while (str_line != ""):
        str_value = str_line.split("(")[1][:-2]
        split_value = str_value.split(",")
        real = float(split_value[0])
        imag = float(split_value[1])
        if (real > pd.max_real): pd.max_real = real
        #    if (config['make_plots'] and real != 0.0):
        #        mag_value = int(math.log10(abs(real)).floor())
        #        if (mag_value > pd.max_real_mag): pd.max_real_mag = mag_value
        #        if (mag_value < pd.min_real_mag): pd.min_real_mag = mag_value
        #        pd.real_hist[mag_value + 99] += 1
        if (imag > pd.max_imag): pd.max_imag = imag
        if (real < pd.min_real): pd.min_real = real
        if (imag < pd.min_imag): pd.min_imag = imag
        cvalue = complex(real, imag)
        absl_val = abs(cvalue)
        if (absl_val > pd.max_absl): pd.max_absl = absl_val
        if (absl_val < pd.sml_absl): pd.sml_absl = absl_val if absl_val > 0.0 else pd.sml_absl
        if (real < 0): real *= -1.0
        if (imag < 0): imag *= -1.0
        if (real < pd.sml_real): pd.sml_real = real if real > 0.0 else pd.sml_real
        if (imag < pd.sml_imag): pd.sml_imag = imag if imag > 0.0 else pd.sml_real
        if (config['make_plots']):
            if (real > 0.0):
                mag_value = int(math.floor(math.log10(real)))
                if (mag_value > pd.max_real_mag): pd.max_real_mag = mag_value
                if (mag_value < pd.min_real_mag): pd.min_real_mag = mag_value
                pd.real_hist[mag_value + 99] += 1
            if (imag > 0.0):
                mag_value = int(math.floor(math.log10(imag)))
                if (mag_value > pd.max_imag_mag): pd.max_imag_mag = mag_value
                if (mag_value < pd.min_imag_mag): pd.min_imag_mag = mag_value
                pd.imag_hist[mag_value + 99] += 1
            if (absl_val > 0.0):
                mag_value = int(math.floor(math.log10(absl_val)))
                if (mag_value > pd.max_absl_mag): pd.max_absl_mag = mag_value
                if (mag_value < pd.min_absl_mag): pd.min_absl_mag = mag_value
                pd.absl_hist[mag_value + 99] += 1
        str_line = matrix_file.readline()
        if (config['limit'] > 0):
            num_read_lines += 1
            if (num_read_lines >= config['limit']):
                str_line = "" # Close the while loop
    matrix_file.close()
    print(" done!")
    pd.log_dynamic_range()
    if (config['make_plots']):
        pd.plot_dynamic_range(config)
    return result

## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of mandatory and optional
#  arguments. The only mandatory argument is the name of the log file to be
#  parsed. Additional optional arguments are an operation filter, which should
#  be the starting sequence of the log strings to pe included in the timing
#  calculation and the number of threads used during code execution.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'help_mode': False,
        'limit': -1,
        'log_plot_y': True,
        'make_plots': False,
        'matrix_name': "",
    }
    for arg in argv[1:]:
        if (arg.startswith("--help")):
            config['help_mode'] = True
        elif (arg.startswith("--limit=")):
            split_arg = arg.split('=')
            config['limit'] = int(split_arg[1])
        elif (arg.startswith("--lin-y")):
            config['log_plot_y'] = False
        elif (arg.startswith("--make-plots")):
            config['make_plots'] = True
        else:
            if (os.path.isfile(arg)):
                config['matrix_name'] = arg
            else:
                raise Exception("Unrecognized argument \'{0:s}\'".format(arg))
    return config

## \brief Print a command-line help summary.
def print_help():
    print("                                                                        ")
    print("****************************   PYDYNRANGE   ****************************")
    print("                                                                        ")
    print("Get the dynamic range of a complex matrix.                              ")
    print("                                                                        ")
    print("Usage: \"./pydynrange.py FILE_NAME [OPTIONS]\"                          ")
    print("                                                                        ")
    print("Valid options are:                                                      ")
    print("--help                                         Print this help and exit.")
    print("--limit=NUMBER                             Check only NUMBER file lines.")
    print("--lin-y                                 Use linear scale on plot y-axes.")
    print("--make-plots        Plot histograms of magnitudes (requires matplotlib).")
    print("                                                                        ")


# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
exit(res)
