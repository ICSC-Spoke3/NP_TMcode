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

## @package pycompare
#  \brief Script to calculate the execution time of logged operations
#
#  Comparing the numeric output can be rendered hard by the amount of information
#  contained in a typical output file and the necessity to determine whether a
#  difference is actually significant or just caused by numeric noise hitting
#  negligible values. The task of `pycompare.py` is to compare two output files, in
#  the assumption that they were written by the FORTRAN and the C++ versions of
#  the code and to flag all the possible inconsistencies according to various
#  severity levels (namely: NOISE, WARNING, and ERROR).
#
#  After execution, the script returns an exit code, which is set to 0, if no
#  error-level inconsistencies were found, or 1 otherwise. This can be used by
#  subsequent system calls to set up a testing suite checking whether the code
#  is able to reproduce legacy results.
#
#  The script execution requires python3.

import re

from sys import argv

## \cond
time_reg = re.compile(r'[0-9]+\.[0-9]+s')
## \endcond

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
        if config['log_name'] is None:
            exit_code = 1
        else:
            operation_time = get_time_from_log(config)
            print("Calculation took %fs."%operation_time)
    return exit_code

## \brief Parse a log file and extract the time.
#
#  \param config: `dict` A dictionary containing the script configuration.
#
#  \returns operation_time: `float` The time of the requested operation in seconds.
def get_time_from_log(config):
    op_time = 0.0
    log_file = open(config['log_name'], 'r')
    file_lines = log_file.readlines()
    log_file.close()
    for li in range(len(file_lines)):
        str_line = file_lines[li]
        if (config['filter'] == "" or str_line.startswith(config['filter'])):
            time_iters = time_reg.finditer(str_line)
            time_groups = []
            for ti in time_iters:
                time_groups.append(ti.group())
            if len(time_groups) == 1:
                op_time += float(time_groups[0][:-1])
    if config['threads'] > 1:
        op_time /= config['threads']
    return op_time

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
        'log_name': None,
        'help_mode': False,
        'filter': "",
        'threads': 1,
    }
    for arg in argv[1:]:
        split_arg = arg.split("=")
        if (arg.startswith("--logname")):
            config['log_name'] = split_arg[1]
        elif (arg.startswith("--filter")):
            config['filter'] = split_arg[1]
        elif (arg.startswith("--threads")):
            config['threads'] = int(split_arg[1])
        elif (arg.startswith("--help")):
            config['help_mode'] = True
        else:
            raise Exception("Unrecognized argument \'{0:s}\'".format(arg))
    return config

## \brief Print a command-line help summary.
def print_help():
    print("                                            ")
    print("***              PYTIMING                ***")
    print("                                            ")
    print("Get the amount of time spent in calculation.")
    print("                                            ")
    print("Usage: \"./pytiming.py OPTIONS\"            ")
    print("                                            ")
    print("Valid options are:                          ")
    print("--logname=TIMING_LOG      File containing log of timing (mandatory).")
    print("--filter=FILTER           Start of the log lines to be accounted for (optional).")
    print("--help                    Print this help and exit.")
    print("--threads=NUM_THREADS     Number of threads or processes used in calculation (optional).")
    print("                                            ")


# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
exit(res)
