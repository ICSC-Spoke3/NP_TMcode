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
#  \brief Script to perform output consistency tests.
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
import os

from math import log10
from sys import argv

## \cond
number_reg = re.compile(r'-?[0-9]\.[0-9]+E?[-+][0-9]{2,5}')
## \endcond

## \brief Main execution code
#
# `main()` is the function that handles the creation of the script configuration
# and the execution of the comparison. It returns an integer value corresponding
# to the number of detected error-level inconsistencies.
#
# \returns errors: `int` Number of detected error-level inconsistencies.
def main():
    config = {}
    try:
        config = parse_arguments()
    except ValueError as ex:
        print(ex)
        print("\nType \"pycompare.py --help\" to get more detailed help.")
        exit(1)
    errors, warnings, noisy = (0, 0, 0)
    if config['help_mode'] or len(argv) == 1:
        config['help_mode'] = True
        print_help()
    else:
        compare_log = compare_files(config)
        errors = compare_log['errors']
        warnings = compare_log['warnings']
        noisy = compare_log['noisy']
        print("ERROR COUNT: %d"%errors)
        print("WARNING COUNT: %d"%warnings)
        print("NOISE COUNT: %d"%noisy)
        if (config['log_html']): reformat_log(config, errors, warnings, noisy)
    if (errors > 0):
        print("FAILURE: {0:s} is not consistent with {1:s}".format(
            config['c_file_name'], config['fortran_file_name']
        ))
    else:
        if (not config['help_mode']):
            print("SUCCESS: {0:s} is consistent with {1:s}".format(
                config['c_file_name'], config['fortran_file_name']
            ))
    return errors

## \brief Perform the comparison of two files.
#
#  The comparison is executed as a line-by-line process. In order to process
#  files correctly, it is required that the two input files have exactly the
#  same format (with the exception of some number formatting subtleties that
#  are handled by regular expressions). Therefore, the first comparison step
#  is testing whether the input files have the same number of lines. If this
#  condition is not met, a formatting problem is assumed and every line in
#  the C++ output file is counted as an error. Otherwise, all the numeric
#  values found in the result are compared for consistency within warning
#  tolerance. Warnings and errors are issued if two values do not match by a
#  fractional amount being, respectively, below or above the threshold for
#  warning. A special case is the detection of suspect numeric noise. This
#  arises on very small quantities as a consequence of differences in the
#  level of approximation or in the hardware implementation of the numeric
#  values, which typically has negligible impact on the overall reslults,
#  even though it can have potentially large fractional mismatches, because
#  it is caused by values that are close to 0.
#
#  Numeric noise is filtered by taking advantage from the fact that the
#  output files are formatted in such a way that values with similar physical
#  meaning are written to the same output line. If the comparison results in
#  a large fractional mismatch on a value that is more than 5 orders of
#  magnitude smaller than the highest order of magnitude that was read from
#  the current row, the discrepancy is flagged as a potential noise effect.
#
#  \param config: `dict` A dictionary containing the script configuration.
#
#  \returns mismatch_count: `tuple(int, int, int)` A tuple that bundles
#  together the numbers of detected errors, warnings and noisy values.
def compare_files(config):
    mismatch_count = {
        'errors': 0,
        'warnings': 0,
        'noisy': 0,
    }
    fortran_file = open(config['fortran_file_name'], 'r')
    c_file = open(config['c_file_name'], 'r')
    l_file = None
    f_lines = []
    c_lines = []
    line_count = 0
    if (not config['linewise']):
        f_lines = fortran_file.readlines()
        c_lines = c_file.readlines()
        line_count = len(f_lines)
        fortran_file.close()
        c_file.close()
    else: # line-wise comparison mode
        f_lines = [fortran_file.readline()]
        c_lines = [c_file.readline()]
        print("INFO: using line-wise mode")
        print("INFO: counting result lines...")
        while (f_lines[0] != ''):
            if (f_lines[0].startswith("INSERTION:")):
                f_lines[0] = fortran_file.readline()
                continue
            if (c_lines[0].startswith("INSERTION:")):
                c_lines[0] = c_file.readline()
            if (c_lines[0] != ''):
                line_count += 1
            else:
                print("ERROR: C++ file is shorter than FORTRAN file.")
                fortran_file.close()
                c_file.close()
                mismatch_count['errors'] = line_count if line_count > 0 else 1
                return mismatch_count
            f_lines[0] = fortran_file.readline()
            c_lines[0] = c_file.readline()
        if (c_lines[0] not in ['', '\n']):
            print("ERROR: C++ file is longer than FORTRAN file.")
            fortran_file.close()
            c_file.close()
            mismatch_count['errors'] = line_count
            return mismatch_count
        fortran_file.close()
        c_file.close()
        print("INFO: the output files have %d lines"%line_count)
        fortran_file = open(config['fortran_file_name'], 'r')
        c_file = open(config['c_file_name'], 'r')
    num_read_lines = 0
    num_filtered_lines = 0
    last_progress = 0;
    # LOG FILE INITIALIZATION #
    if (config['log_html']):
        l_file = open(config['html_output'], 'w')
        l_file.write("<!DOCTYPE html>\n")
        l_file.write("<html xmnls=\"http://www.w3.org/1999/xhtml\">\n")
        l_file.write("  <header>\n")
        l_file.write(
            "    <h1>Comparison between {0:s} and {1:s}</h1>\n".format(
                config['fortran_file_name'], config['c_file_name']
            )
        )
        l_file.write("  </header>\n")
        l_file.write("  <body>\n")
        l_file.write("    <div>Numeric noise is marked <span style=\"font-weight: bold; color: rgb(0,185,0)\">"
                     + "GREEN</span>, warnings are marked <span style=\"font-weight: bold; color: rgb(0,0,255)\">"
                     + "BLUE</span> and errors are marked <span style=\"font-weight: bold; color: rgb(255,0,0)\">"
                     + "RED</span>.</div>\n")
    # END LOG FILE INITIALIZATION #
    line_loop = True
    num_len = 1
    if (line_count > 0):
        num_len = max(4, int(log10(line_count)) + 1)
    if (config['say_progress']):
        print("INFO: checking file contents...   0", end='%', flush=True)
    else:
        print("INFO: checking file contents...")
    while (line_loop):
        if (not config['linewise']):
            line_loop = False
        else:
            f_lines[0] = fortran_file.readline()
            if (f_lines[0].startswith("INSERTION:")):
                continue
            c_lines[0] = c_file.readline()
            if (c_lines[0].startswith("INSERTION:")):
                c_lines = [c_file.readline()]
                num_read_lines += 1
            num_read_lines += 1
            num_filtered_lines += 1
        # Start here the comparison loop
        if (len(f_lines) == len(c_lines)):
            for li in range(len(f_lines)):
                line_result = compare_lines(f_lines[li], c_lines[li], config, num_read_lines, num_len, l_file)
                mismatch_count['errors'] += line_result[0]
                mismatch_count['warnings'] += line_result[1]
                mismatch_count['noisy'] += line_result[2]
                if (mismatch_count['errors'] > 0 and not config['check_all']):
                    print("\nINFO: mismatch found at line %d"%(num_read_lines))
                    line_loop = False
                    break
        else:
            mismatch_count['errors'] = len(c_lines)
            print("\nERROR: {0:s} and {1:s} have different numbers of lines!".format(
                config['fortran_file_name'], config['c_file_name']
            ))
            if (config['log_html']):
                print("Different file sizes. No log produced.")
                config['log_html'] = False
        if (num_filtered_lines >= line_count):
            line_loop = False
        if (config['say_progress']):
            progress = int(100 * num_filtered_lines / line_count)
            if (progress > last_progress):
                print("\b\b\b\b%3d"%progress, end="%", flush=True)
                last_progress = progress
        #End line loop
    if l_file is not None:
        l_file.write("  </body>\n")
        l_file.write("</html>\n")
        l_file.close()
    if (config['say_progress']):
        print("")
    return mismatch_count

## \brief Perform the comparison of two file lines.
#
#  This function handles the line-by-line comparison of coded result files. Depending
#  on whether a HTML log report was requested, it also undertakes the task of
#  formatting the HTML code to show the comparison results as high-lighted entries,
#  according to the severity degree of the mismatch.
#
#  \param f_line: `string` A line extracted from the FORTRAN output file.
#  \param c_line: `string` A line extracted from the C++ output file.
#  \param config: `dict` A dictionary containing the script configuration.
#  \param line_num: `int` The number of the current line (0-indexed).
#  \param num_len: `int` The number digits to format the line number tag in the HTML log.
#  \param log_file: `file` A file where to write logging information, if required.
#
#  \returns mismatch_count: `tuple(int, int, int)` A tuple that bundles
#  together the numbers of detected errors, warnings and noisy values.
def compare_lines(f_line, c_line, config, line_num=0, num_len=4, log_file=None):
    errors = 0
    warnings = 0
    noisy = 0
    f_line = f_line.replace("D-","E-").replace("D+","E+")
    ref_format = "    <div><span style=\"font-weight: bold; color: rgb(125,125,125)\"><pre><code>{0:%ds}"%num_len
    ref_line = (ref_format + ": {1:s}</code></pre></span></div>\n").format("ORIG", f_line[:-1])
    log_line = ""
    if (f_line == c_line):
        if log_file is not None:
            if (config['full_log']):
                num_format = "    <div><pre><code>{0:0%dd}"%num_len
                log_line = (num_format + ": {1:s}</code></pre></div>\n").format(line_num, c_line[:-1])
                log_file.write(log_line)
    else:
        iter_f_values = number_reg.finditer(f_line)
        iter_c_values = number_reg.finditer(c_line)
        f_starts, f_ends, f_groups = [], [], []
        c_starts, c_ends, c_groups = [], [], []
        for fi in iter_f_values:
            f_starts.append(fi.start())
            f_ends.append(fi.end())
            f_groups.append(fi.group())
        for ci in iter_c_values:
            c_starts.append(ci.start())
            c_ends.append(ci.end())
            c_groups.append(ci.group())

        if (len(f_groups) == len(c_groups)):
            severities = mismatch_severities(f_groups, c_groups, config)
            if log_file is not None:
                if (len(severities) > 0):
                    num_format = "    <div><pre><code>{0:0%dd}"%num_len
                    log_line = (num_format + ": ").format(line_num)
                    log_line = log_line + c_line[0:c_starts[0]]
            for si in range(len(severities) - 1):
                if (severities[si] == 1): noisy += 1
                elif (severities[si] == 2): warnings += 1
                elif (severities[si] == 3): errors += 1
                if log_file is not None:
                    if (severities[si] == 0):
                        log_line = log_line + c_groups[si] + c_line[c_ends[si]:c_starts[si + 1]]
                    elif (severities[si] == 1):
                        log_line = (
                            log_line + "</code><span style=\"font-weight: bold; color: rgb(0,185,0)\"><code>"
                            + c_groups[si] + "</code></span><code>" + c_line[c_ends[si]:c_starts[si + 1]]
                        )
                    elif (severities[si] == 2):
                        log_line = (
                            log_line + "</code><span style=\"font-weight: bold; color: rgb(0,0,255)\"><code>"
                            + c_groups[si] + "</code></span><code>" + c_line[c_ends[si]:c_starts[si + 1]]
                        )
                    elif (severities[si] == 3):
                        log_line = (
                            log_line + "</code><span style=\"font-weight: bold; color: rgb(255,0,0)\"><code>"
                            + c_groups[si] + "</code></span><code>" + c_line[c_ends[si]:c_starts[si + 1]]
                        )
            if (len(severities) > 0):
                # Single error test modification
                if (severities[-1] == 1): noisy += 1
                elif (severities[-1] == 2): warnings += 1
                elif (severities[-1] == 3):
                    split_c_line = c_line.split('/')
                    if (config['warning_threshold'] == 0.0): errors += 1
                    elif (len(split_c_line) != 2): errors += 1
            if log_file is not None:
                if (len(severities) > 0):
                    if (severities[-1] == 0):
                        log_line = (
                            log_line + c_groups[-1] + c_line[c_ends[-1]:len(c_line) - 1]
                        )
                    if (severities[-1] == 1):
                        log_line = (
                            log_line + "</code><span style=\"font-weight: bold; color: rgb(0,185,0)\"><code>"
                            + c_groups[-1] + "</code></span><code>" + c_line[c_ends[-1]:len(c_line) - 2]
                        )
                    if (severities[-1] == 2):
                        log_line = (
                            log_line + "</code><span style=\"font-weight: bold; color: rgb(0,0,255)\"><code>"
                            + c_groups[-1] + "</code></span><code>" + c_line[c_ends[-1]:len(c_line) - 2]
                        )
                    if (severities[-1] == 3):
                        split_c_line = c_line.split('/')
                        if (len(split_c_line) == 2):
                            if (config['warning_threshold'] != 0.0):
                                log_line = (
                                    log_line + "</code><span style=\"font-weight: bold; color: rgb(0,185,0)\"><code>"
                                    + c_groups[-1] + "</code></span><code>" + c_line[c_ends[-1]:len(c_line) - 2]
                                )
                            else:
                                log_line = (
                                    log_line + "</code><span style=\"font-weight: bold; color: rgb(255,0,0)\"><code>"
                                    + c_groups[-1] + "</code></span><code>" + c_line[c_ends[-1]:len(c_line) - 2]
                                )
                        else:
                            log_line = (
                                log_line + "</code><span style=\"font-weight: bold; color: rgb(255,0,0)\"><code>"
                                + c_groups[-1] + "</code></span><code>" + c_line[c_ends[-1]:len(c_line) - 2]
                            )
                    if ((not config['hide_noise'] and noisy > 0) or warnings > 0 or errors > 0):
                        log_file.write(log_line + "</code></pre></div>\n")
                        log_file.write(ref_line)
        else: # The two lines contain a different number of numeric values
            if (log_file is not None):
                num_format = "    <div><pre><code>{0:0%dd}"%num_len
                log_line = (num_format + ": ").format(line_num)
                log_line = (
                    log_line + "</code><span style=\"font-weight: bold; color: rgb(255,0,0)\"><code>"
                    + c_line + "</code></span><code>"
                )
                log_file.write(log_line + "</code></pre></div>\n")
                log_file.write(ref_line)
            errors += 1
    return (errors, warnings, noisy)

## \brief Determine the severity of a numerical mismatch.
#
#  The severity scale is currently designed with the following integer codes:
#
#  0 - the values are equal
#
#  1 - the values are subject to suspect numerical noise (green fonts)
#
#  2 - the values are different but below error threshold (blue fonts)
#
#  3 - the values differ more than error threshold (red fonts)
#
#  \param str_f_values: `array(string)` The strings representing the numeric
#     values read from the FORTRAN output file.
#  \param str_c_values: `array(string)` The strings representing the numeric
#     values read from the C++ output file.
#  \param config: `dict` A dictionary containing the configuration options from
#     which to read the warning and the error threshold.
#
#  \returns result: `array(int)` An array of severity codes ordered as the
#  input numeric values.
def mismatch_severities(str_f_values, str_c_values, config):
    result = []
    if len(str_f_values) == len(str_c_values):
        result = [0 for ri in range(len(str_f_values))]
        f_values = []
        c_values = []
        # Convert numeric strings to numbers
        for i in range(len(str_f_values)):
            # Add the exponent marker if it is missing
            temp_str_value = str_f_values[i][1:]
            split_temp = temp_str_value.split('-')
            if len(split_temp) > 1:
                if (split_temp[0][-1] != 'E'):
                    str_f_values[i] = str_f_values[i][0] + split_temp[0] + "E-" + split_temp[1]
            f_values.append(float(str_f_values[i]))
            temp_str_value = str_c_values[i][1:]
            split_temp = temp_str_value.split('-')
            if len(split_temp) > 1:
                if (split_temp[0][-1] != 'E'):
                    str_c_values[i] = str_c_values[i][0] + split_temp[0] + "E-" + split_temp[1]
            c_values.append(float(str_c_values[i]))
            # End of missing exponent marker correction
        # End string to number conversion
        # Evaluate the maximum scale
        max_f_log = -1.0e12
        max_c_log = -1.0e12
        for si in range(len(f_values)):
            if (f_values[si] != 0):
                sign = 1.0 if f_values[si] > 0.0 else -1.0
                log_f_value = log10(sign * f_values[si])
                if (log_f_value > max_f_log): max_f_log = log_f_value
            if (c_values[si] != 0):
                sign = 1.0 if c_values[si] > 0.0 else -1.0
                log_c_value = log10(sign * c_values[si])
                if (log_c_value > max_c_log): max_c_log = log_c_value
        if (max_f_log == -1.0e12): max_f_log = 0.0
        if (max_c_log == -1.0e12): max_c_log = 0.0
        # End of maximum scale evaluation
        # Compare the numbers
        for i in range(len(f_values)):
            if (f_values[i] != c_values[i]):
                if (f_values[i] != 0.0):
                    sign = 1.0 if f_values[i] > 0.0 else -1.0
                    log_f_value = log10(sign * f_values[i])
                    if (log_f_value > max_f_log - 5.0):
                        scale = 10.0**(log_f_value - max_f_log)
                        fractional = scale * (f_values[i] - c_values[i]) / f_values[i]
                        if (fractional < 0.0): fractional *= -1.0
                        if (fractional <= config['warning_threshold']):
                            if (log_f_value > config['data_order']):
                                result[i] = 2
                            else:
                                result[i] = 1
                        else:
                            if (log_f_value > config['data_order']):
                                result[i] = 3
                            else:
                                result[i] = 1
                    else:
                        result[i] = 1
                else: # f_values[i] == 0 and c_values[i] != 0
                    sign = 1.0 if c_values[i] > 0.0 else -1.0
                    log_c_value = log10(sign * c_values[i])
                    if (log_c_value > max_c_log - 5.0):
                        scale = 10.0**(log_c_value - max_c_log)
                        fractional = scale * (c_values[i] - f_values[i]) / c_values[i]
                        if (fractional < 0.0): fractional *= -1.0
                        if (fractional <= config['warning_threshold']):
                            if (log_c_value > config['data_order']):
                                result[i] = 2
                            else:
                                result[i] = 1
                        else:
                            if (log_c_value > config['data_order']):
                                result[i] = 3
                            else:
                                result[i] = 1
                    else:
                        result[i] = 1
        # End number comparison
    return result
    
## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of mandatory and optional
#  arguments. Mandatory arguments are those required to execute a meaningful
#  comparison and they are limited to the names of the files that need to be
#  compared. The other arguments affect whether the script should produce an
#  HTML log file and what level of detail needs to be included in this log.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'c_file_name': '',
        'check_all': True,
        'data_order': -99.0,
        'fortran_file_name': '',
        'full_log': False,
        'help_mode': False,
        'hide_noise': True,
        'html_output': 'pycompare.html',
        'linewise': True,
        'log_html': False,
        'say_progress': True,
        'warning_threshold': 0.005,
    }
    arg_index = 1
    skip_arg = False
    for arg in argv[1:]:
        if skip_arg:
            skip_arg = False
            continue
        split_arg = arg.split("=")
        if (arg.startswith("--ffile")):
            config['fortran_file_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--cfile")):
            config['c_file_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--data-order")):
            config['data_order'] = float(split_arg[1])
        elif (arg.startswith("--full")):
            config['full_log'] = True
        elif (arg.startswith("--html")):
            config['log_html'] = True
            if (len(split_arg) == 2):
                config['html_output'] = split_arg[1]
        elif (arg.startswith("--help")):
            config['help_mode'] = True
        elif (arg.startswith("--linewise")):
            config['linewise'] = True
        elif (arg.startswith("--no-progress")):
            config['say_progress'] = False
        elif (arg.startswith("--quick")):
            config['check_all'] = False
        elif (arg.startswith("--show-noise")):
            config['hide_noise'] = False
        elif (arg.startswith("--warn")):
            config['warning_threshold'] = float(split_arg[1])
        else:
            raise ValueError("Unrecognized argument \'{0:s}\'".format(arg))
        arg_index += 1
    return config

## \brief Print a command-line help summary.
def print_help():
    print("                                            ")
    print("***              PYCOMPARE               ***")
    print("                                            ")
    print("Compare the output of C++ and FORTRAN codes.")
    print("                                            ")
    print("Usage: \"./pycompare.py OPTIONS\"           ")
    print("                                            ")
    print("Valid options are:                          ")
    print("--ffile FORTRAN_OUTPUT    File containing the output of the FORTRAN code (mandatory).")
    print("--cfile C++_OUTPUT        File containing the output of the C++ code (mandatory).")
    print("--data-order=ORDER        Consider data only down to specified order (default order is -99).")
    print("--full                    Print all lines to log file (default prints only mismatches).")
    print("--help                    Print this help and exit.")
    print("--html[=OPT_OUTPUT_NAME]  Enable logging to HTML file (default logs to \"pycompare.html\").")
    print("--linewise                Load only one line at a time. Useful to compare big files (True by default).")
    print("--no-progress             Disable progress logging.")
    print("--quick                   Stop on first mismatch (default is to perform a full check).")
    print("--show-noise              Show noise in reports (default is to hide it).")
    print("--warn                    Set a fractional threshold for numeric warning (default = 0.005).")
    print("                                            ")

## \brief Add summary information to the HTML log file
#
#  In the case when a HTML log is requested, it is useful to obtain an overview
#  of the detected inconsistencies. This function undertakes the task of adding
#  a summary of the error, warning and noise counts on top of the log.
#  
#  \param config: `dict` A dictionary containing the script configuration.
#  \param errors: `int` The number of errors detected by the comparison.
#  \param warnings: `int` The number of warnings detected by the comparison.
#  \param noisy: `int` The number of noisy values detected by the comparison.
def reformat_log(config, errors, warnings, noisy):
    log_file = open(config['html_output'], 'r')
    new_file = open("PYCOMPARE_TEMPORARY_LOG.html", 'w')
    for hi in range(7):
        log_line = log_file.readline()
        new_file.write(log_line)
    str_errors = "error" if errors == 1 else "errors"
    str_warnings = "warning" if warnings == 1 else "warnings"
    str_noisy = "noisy value" if noisy == 1 else "noisy values"
    summary = "    <div>Comparison yielded %d %s"%(errors, str_errors)
    summary = summary + ", %d %s"%(warnings, str_warnings)
    summary = summary + " and %d %s.</div>\n"%(noisy, str_noisy)
    new_file.write(summary)
    log_line = log_file.readline()
    while (log_line != ''):
        new_file.write(log_line)
        log_line = log_file.readline()
    log_file.close()
    new_file.close()
    os.remove(config['html_output'])
    os.rename("PYCOMPARE_TEMPORARY_LOG.html", config['html_output'])

# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
if (res > 0): exit(1)
exit(0)
