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

## @package parse_output
#  \brief Script to extract physical quantities from NPTMcode output files.
#
#  This script is intended to assist users in the inspection of results. It
#  can be used to extract physical quantities, such as radiation pressure
#  forces and cross-sections from the standard code output files and save
#  the data in machine readable CSV files.
#
#  The script execution requires python3.

import math
from sys import argv

## \brief Main execution code
#
# `main()` is the function that handles the creation of the script configuration
# and the execution of the comparison. It returns an integer value corresponding
# to the number of detected error-level inconsistencies.
#
# \returns errors: `int` Number of runtime errors (0 means successful run).
def main():
    config = {}
    errors = 0
    try:
        config = parse_arguments()
    except ValueError as ex:
        print(ex)
        print("\nType \"parse_output.py --help\" to get more detailed help.")
        errors = 1
    if config['help_mode']:
        config['help_mode'] = True
        print_help()
    else:
        if (config['result_file'] == ''):
            print("ERROR: no input file name given (missing \"--in FILE\" option).") 
            errors += 1
        if (config['output_name'] == ''):
            print("ERROR: no output file name given (missing \"--out FILE\" option).")
            errors += 1
        if (config['application'] == ''):
            print("ERROR: no application identifier given (missing \"--app=APP\" option).")
            errors += 1
        if (errors == 0):
            if (config['format'] == 'LEGACY'):
                if (config['application'] == 'CLU'):
                    try:
                        errors += parse_legacy_oclu(config)
                    except Exception as ex:
                        print(ex)
                        errors += 1
                if (config['application'] == 'SPH'):
                    try:
                        errors += parse_legacy_osph(config)
                    except Exception as ex:
                        print(ex)
                        errors += 1
    return errors

## \brief Parse the command line arguments.
#
#  The script behaviour can be modified through a set of mandatory and optional
#  arguments. Mandatory arguments are those required to execute a meaningful
#  parsing and they are limited to the names of the files that need to be read
#  and written. The other arguments affect the format and the data fields that
#  are sought for.
#
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments():
    config = {
        'result_file': '',
        'output_name': '',
        'application': '',
        'format': 'LEGACY',
        'selection': '[ALL]',
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
            config['result_file'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
            if (config['result_file'].endswith('.hd5')):
                config['format'] = 'HDF5'
        elif (arg.startswith("--out")):
            config['output_name'] = argv[arg_index + 1]
            arg_index += 1
            skip_arg = True
        elif (arg.startswith("--app=")):
            config['application'] = split_arg[1]
        elif (arg.startswith("--selection=")):
            config['selection'] = split_arg[1]
        elif (arg.startswith("--help")):
            config['help_mode'] = True
        else:
            raise ValueError("Unrecognized argument \'{0:s}\'".format(arg))
        arg_index += 1
    return config

## \brief Parse a legacy output file of np_cluster.
#
#  \param config: `dict` A dictionary containing the script configuration.
#  \return errors: `int` The number of encountered errors.
def parse_legacy_oclu(config):
    errors = 0
    oclu_name = config['result_file']
    root_name = config['output_name']
    out20 = None # Cluster average cross-sections
    out31 = None # Cluster differential cross-sections in state -1
    out32 = None # Cluster differential cross-sections in state +1
    out40 = None # Cluster integrated asymmetry parameter and radiation pressure
    out51 = None # Cluster differential asymmetry parameter and radiation pressure forces in state -1
    out52 = None # Cluster differential asymmetry parameter and radiation pressure forces in state +1
    try:
        oclu_file = open(oclu_name, "r") # open the OCLU file for reading
        file_line = "first line" # a string to parse the OCLU file lines
        if ('ALL' in config['selection'] or 'ICS' in config['selection']):
            out20 = open(root_name + "_ics.csv", "w") # open a file for integrated cross sections
            out20.write("Wavelength,ScaSec,AbsSec,ExtSec\n")
        if ('ALL' in config['selection'] or 'DCS' in config['selection']):
            out31 = open(root_name + "_dcs1.csv", "w") # open a file for differential cross-sections in state -1
            out31.write("Wavelength,THi,THs,PHi,PHs,ScaSec,AbsSec,ExtSec\n")
            out32 = open(root_name + "_dcs2.csv", "w") # open a file for differential cross-sections in state +1
            out32.write("Wavelength,THi,THs,PHi,PHs,ScaSec,AbsSec,ExtSec\n")
        if ('ALL' in config['selection'] or 'IRP' in config['selection']):
            out40 = open(root_name + "_irp.csv", "w") # open a file for integrated radiation pressure forces
            out40.write("Wavelength,CosAv,RaPr\n")
        if ('ALL' in config['selection'] or 'DRP' in config['selection']):
            out51 = open(root_name + "_drp1.csv", "w") # open a file for differential radiation pressure forces in state -1
            out51.write("Wavelength,THi,THs,PHi,PHs,CosAv,RaPr,Fl,Fr,Fk,Fx,Fy,Fz,TQEl,TQEr,TQEk,TQEx,TQEy,TQEz,TQSl,TQSr,TQSk,TQSx,TQSy,TQSz\n")
            out52 = open(root_name + "_drp2.csv", "w") # open a file for differential radiation pressure forces in state +1
            out52.write("Wavelength,THi,THs,PHi,PHs,CosAv,RaPr,Fl,Fr,Fk,Fx,Fy,Fz,TQEl,TQEr,TQEk,TQEx,TQEy,TQEz,TQSl,TQSr,TQSk,TQSx,TQSy,TQSz\n")

        # Define the quantities that you need to extract
        alam = 0.0
        vk = 0.0
        scaleOnXi = False

        # Read the output file preamble
        for i in range(2):
            file_line = oclu_file.readline()
        nsph = int(file_line[0:6])
        for i in range(nsph + 3):
            file_line = oclu_file.readline()
        thifirst = float(file_line[0:11].replace("D", "E"))
        thistep = float(file_line[12:21].replace("D", "E"))
        thilast = float(file_line[22:31].replace("D", "E"))
        thsfirst = float(file_line[32:41].replace("D", "E"))
        thsstep = float(file_line[42:51].replace("D", "E"))
        thslast = float(file_line[52:61].replace("D", "E"))
        nthi = 1 if thistep == 0.0 else 1 + int((thilast - thifirst) / thistep)
        nths = 1 if thsstep == 0.0 else 1 + int((thslast - thsfirst) / thsstep)
        for i in range(2):
            file_line = oclu_file.readline()
        phifirst = float(file_line[0:11].replace("D", "E"))
        phistep = float(file_line[12:21].replace("D", "E"))
        philast = float(file_line[22:31].replace("D", "E"))
        phsfirst = float(file_line[32:41].replace("D", "E"))
        phsstep = float(file_line[42:51].replace("D", "E"))
        phslast = float(file_line[52:61].replace("D", "E"))
        nphi = 1 if phistep == 0.0 else 1 + int((philast - phifirst) / phistep)
        nphs = 1 if phsstep == 0.0 else 1 + int((phslast - phsfirst) / phsstep)
        ndirs = nthi * nths * nphi * nphs
        
        while ("JXI =" not in file_line):
            file_line = oclu_file.readline() # read the next OSPH file line
            if ("XI IS SCALE FACTOR FOR LENGTHS" in file_line):
                scaleOnXi = True
                vk = float(file_line[5:20].replace('D', 'E'))
        
        # Parsing loop until the end of the OCLU file
        while (file_line != ""):
            file_line = oclu_file.readline() # read the next OCLU file line
            if (scaleOnXi):
                if (file_line.startswith("  XI=")):
                    xi = float(file_line[5:20].replace('D', 'E'))
                    alam = 2.0 * math.pi * xi / vk
            else:
                if (file_line.startswith("  VK=")):
                    # we found VK, so we calculate lambda
                    # extract VK as a number from a string section, after
                    # replacing FORTRAN's D with E
                    vk = float(file_line[5:20].replace("D", "E"))
                    alam = 2.0 * math.pi / vk
            if ("CLUSTER (ENSEMBLE AVERAGE, MODE 0)" in file_line):
                # we are in average section. We start a nested loop to
                # extract the average values
                found_averages = False
                while (not found_averages):
                    file_line = oclu_file.readline()
                    if ("----- SCC ----- ABC ----- EXC ----- ALBEDC --" in file_line):
                        # we know we are in LIN -1 because it is the first one
                        # we also know that the next line contains the values
                        # we are looking for, so we read it
                        file_line = oclu_file.readline()
                        # we now extract the values from string sections
                        scasm = float(file_line[1:15].replace("D", "E"))
                        abssm = float(file_line[17:30].replace("D", "E"))
                        extsm = float(file_line[32:45].replace("D", "E"))
                        # we can write the average values, similarly to fort.22
                        # note that \n puts a new line at the end of the string
                        output_line = "{0:.7E},{1:.7E},{2:.7E},{3:.7E}\n".format(alam, scasm, abssm, extsm)
                        if (out20 is not None): out20.write(output_line)
                        # we know that the asymmetry parameter and the radiation
                        # pressure forces will be after 8 more lines
                        for i in range(8):
                            file_line = oclu_file.readline()
                        cosav = float(file_line[8:23].replace("D", "E"))
                        rapr = float(file_line[31:46].replace("D", "E"))
                        output_line = "{0:.7E},{1:.7E},{2:.7E}\n".format(alam, cosav, rapr)
                        if (out40 is not None): out40.write(output_line)
                        found_averages = True # terminate the inner loop
                # the averages were written. We look for CLUSTER section
                # using another inner loop
                found_differentials = False
                while (not found_differentials):
                    for di in range(ndirs):
                        while ("JTH =" not in file_line):
                            file_line = oclu_file.readline()
                        file_line = oclu_file.readline()
                        tidg = float(file_line[7:17].replace("D", "E"))
                        pidg = float(file_line[24:34].replace("D", "E"))
                        tsdg = float(file_line[41:51].replace("D", "E"))
                        psdg = float(file_line[58:68].replace("D", "E"))
                        while ("     CLUSTER" not in file_line):
                            file_line = oclu_file.readline()
                        if ("     CLUSTER" in file_line):
                            # we found CLUSTER. We know cross-sections for
                            # polarization state -1 will be after 2 more lines
                            for i in range(3):
                                file_line = oclu_file.readline()
                                # the following check is needed to parse C++ output
                                if ("INSERTION" in file_line):
                                    file_line = oclu_file.readline()
                            scc31 = float(file_line[1:15].replace("D", "E"))
                            abc31 = float(file_line[17:30].replace("D", "E"))
                            exc31 = float(file_line[32:45].replace("D", "E"))
                            # we can write the differential values, similarly to fort.31
                            output_line = "{0:.7E},{1:.3E},{2:.3E},{3:.3E},{4:.3E},{5:.7E},{6:.7E},{7:.7E}\n".format(alam, tidg, tsdg, pidg, psdg, scc31, abc31, exc31)
                            if (out31 is not None): out31.write(output_line)
                            # we know that RAPRS values for polarization state -1
                            # are after 9 more lines
                            for i in range(9):
                                file_line = oclu_file.readline()
                                # the following check is needed to parse C++ output
                                if ("INSERTION" in file_line):
                                    file_line = oclu_file.readline()
                            cosav = float(file_line[8:23].replace("D", "E"))
                            rapr = float(file_line[31:46].replace("D", "E"))
                            # we read the forces and torques
                            file_line = oclu_file.readline()
                            fl = float(file_line[5:20].replace("D", "E"))
                            fr = float(file_line[25:40].replace("D", "E"))
                            fk = float(file_line[45:60].replace("D", "E"))
                            file_line = oclu_file.readline()
                            fx = float(file_line[5:20].replace("D", "E"))
                            fy = float(file_line[25:40].replace("D", "E"))
                            fz = float(file_line[45:60].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQEl = float(file_line[8:23].replace("D", "E"))
                            TQEr = float(file_line[31:46].replace("D", "E"))
                            TQEk = float(file_line[54:69].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQSl = float(file_line[8:23].replace("D", "E"))
                            TQSr = float(file_line[31:46].replace("D", "E"))
                            TQSk = float(file_line[54:69].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQEx = float(file_line[8:23].replace("D", "E"))
                            TQEy = float(file_line[31:46].replace("D", "E"))
                            TQEz = float(file_line[54:69].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQSx = float(file_line[8:23].replace("D", "E"))
                            TQSy = float(file_line[31:46].replace("D", "E"))
                            TQSz = float(file_line[54:69].replace("D", "E"))
                            # we can write the RAPRS values
                            output_line = "{0:.7E},{1:.3E},{2:.3E},{3:.3E},{4:.3E},{5:.7E},{6:.7E},{7:.7E},{8:.7E},{9:.7E},{10:.7E},{11:.7E},{12:.7E},{13:.7E},{14:.7E},{15:.7E},{16:.7E},{17:.7E},{18:.7E},{19:.7E},{20:.7E},{21:.7E},{22:.7E},{23:.7E},{24:.7E}\n".format(alam, tidg, tsdg, pidg, psdg, cosav, rapr, fl, fr, fk, fx, fy, fx, TQEl, TQEr, TQEk, TQEx, TQEy, TQEz, TQSl, TQSr, TQSk, TQSx, TQSy, TQSz)
                            if (out51 is not None): out51.write(output_line)
                            # we know the differential values for polarization
                            # state 1 are after 3 more lines
                            for i in range(3):
                                file_line = oclu_file.readline()
                                # the following check is needed to parse C++ output
                                if ("INSERTION" in file_line):
                                    file_line = oclu_file.readline()
                            scc32 = float(file_line[1:15].replace("D", "E"))
                            abc32 = float(file_line[17:30].replace("D", "E"))
                            exc32 = float(file_line[32:45].replace("D", "E"))
                            # we can write the differential values, similarly to fort.31
                            output_line = "{0:.7E},{1:.3E},{2:.3E},{3:.3E},{4:.3E},{5:.7E},{6:.7E},{7:.7E}\n".format(alam, tidg, tsdg, pidg, psdg, scc32, abc32, exc32)
                            if (out32 is not None): out32.write(output_line)
                            # we know that RAPRS values for polarization state 1
                            # are after 9 more lines
                            for i in range(9):
                                file_line = oclu_file.readline()
                                # the following check is needed to parse C++ output
                                if ("INSERTION" in file_line):
                                    file_line = oclu_file.readline()
                            cosav = float(file_line[8:23].replace("D", "E"))
                            rapr = float(file_line[31:46].replace("D", "E"))
                            # we read the forces and torques
                            file_line = oclu_file.readline()
                            fl = float(file_line[5:20].replace("D", "E"))
                            fr = float(file_line[25:40].replace("D", "E"))
                            fk = float(file_line[45:60].replace("D", "E"))
                            file_line = oclu_file.readline()
                            fx = float(file_line[5:20].replace("D", "E"))
                            fy = float(file_line[25:40].replace("D", "E"))
                            fz = float(file_line[45:60].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQEl = float(file_line[8:23].replace("D", "E"))
                            TQEr = float(file_line[31:46].replace("D", "E"))
                            TQEk = float(file_line[54:69].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQSl = float(file_line[8:23].replace("D", "E"))
                            TQSr = float(file_line[31:46].replace("D", "E"))
                            TQSk = float(file_line[54:69].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQEx = float(file_line[8:23].replace("D", "E"))
                            TQEy = float(file_line[31:46].replace("D", "E"))
                            TQEz = float(file_line[54:69].replace("D", "E"))
                            file_line = oclu_file.readline()
                            TQSx = float(file_line[8:23].replace("D", "E"))
                            TQSy = float(file_line[31:46].replace("D", "E"))
                            TQSz = float(file_line[54:69].replace("D", "E"))
                            # we can write the RAPRS values
                            output_line = "{0:.7E},{1:.3E},{2:.3E},{3:.3E},{4:.3E},{5:.7E},{6:.7E},{7:.7E},{8:.7E},{9:.7E},{10:.7E},{11:.7E},{12:.7E},{13:.7E},{14:.7E},{15:.7E},{16:.7E},{17:.7E},{18:.7E},{19:.7E},{20:.7E},{21:.7E},{22:.7E},{23:.7E},{24:.7E}\n".format(alam, tidg, tsdg, pidg, psdg, cosav, rapr, fl, fr, fk, fx, fy, fx, TQEl, TQEr, TQEk, TQEx, TQEy, TQEz, TQSl, TQSr, TQSk, TQSx, TQSy, TQSz)
                            if (out52 is not None): out52.write(output_line)
                    found_differentials = True # terminate the inner loop
            # The parsing loop ends here

        if (out20 is not None): out20.close()
        if (out31 is not None): out31.close()
        if (out32 is not None): out32.close()
        if (out40 is not None): out40.close()
        if (out51 is not None): out51.close()
        if (out52 is not None): out52.close()
        oclu_file.close() # close the OCLU file
    except Exception as ex:
        print(ex)
        errors += 1
    return errors

## \brief Parse a legacy output file of np_sphere.
#
#  \param config: `dict` A dictionary containing the script configuration.
#  \return errors: `int` The number of encountered errors.
def parse_legacy_osph(config):
    errors = 0
    osph_name = config['result_file']
    root_name = config['output_name']
    out20 = None # Sphere average cross-sections
    out30 = None # Sphere integrated asymmetry parameter and radiation pressure
    out40 = None # Sphere differential radiation forces
    try:
        osph_file = open(osph_name, "r") # open the OSPH file for reading
        file_line = "first line" # a string to parse the OSPH file lines
        if ('ALL' in config['selection'] or 'ICS' in config['selection']):
            out20 = open(root_name + "_ics.csv", "w") # open a file for integrated cross sections
            out20.write("Wavelength,ScaSec,AbsSec,ExtSec\n")
        if ('ALL' in config['selection'] or 'IRP' in config['selection']):
            out30 = open(root_name + "_irp.csv", "w") # open a file for integrated radiation pressure forces
            out30.write("Wavelength,CosAv,RaPr\n")
        if ('ALL' in config['selection'] or 'DRP' in config['selection']):
            out40 = open(root_name + "_drp.csv", "w") # open a file for differential radiation pressure forces in state -1
            out40.write("Wavelength,THi,THs,PHi,PHs,Fx,Fy,Fz\n")

        # Define the quantities that you need to extract
        alam = 0.0
        scaleOnXi = False
        vk = 0.0

        # Read the output file preamble
        for i in range(2):
            file_line = osph_file.readline()
        nsph = int(file_line[0:6])
        if (nsph != 1):
            errors += 1
            print("ERROR: number of spheres is not 1!")
            return errors
        for i in range(2):
            file_line = osph_file.readline()
        thifirst = float(file_line[0:11].replace("D", "E"))
        thistep = float(file_line[12:21].replace("D", "E"))
        thilast = float(file_line[22:31].replace("D", "E"))
        thsfirst = float(file_line[32:41].replace("D", "E"))
        thsstep = float(file_line[42:51].replace("D", "E"))
        thslast = float(file_line[52:61].replace("D", "E"))
        nthi = 1 if thistep == 0.0 else 1 + int((thilast - thifirst) / thistep)
        nths = 1 if thsstep == 0.0 else 1 + int((thslast - thsfirst) / thsstep)
        for i in range(2):
            file_line = osph_file.readline()
        phifirst = float(file_line[0:11].replace("D", "E"))
        phistep = float(file_line[12:21].replace("D", "E"))
        philast = float(file_line[22:31].replace("D", "E"))
        phsfirst = float(file_line[32:41].replace("D", "E"))
        phsstep = float(file_line[42:51].replace("D", "E"))
        phslast = float(file_line[52:61].replace("D", "E"))
        nphi = 1 if phistep == 0.0 else 1 + int((philast - phifirst) / phistep)
        nphs = 1 if phsstep == 0.0 else 1 + int((phslast - phsfirst) / phsstep)
        ndirs = nthi * nths * nphi * nphs

        while ("JXI =" not in file_line):
            file_line = osph_file.readline() # read the next OSPH file line
            if ("XI IS SCALE FACTOR FOR LENGTHS" in file_line):
                scaleOnXi = True
                vk = float(file_line[5:20].replace('D', 'E'))
        
        # Parsing loop until the end of the OSPH file
        while (file_line != ""):
            file_line = osph_file.readline() # read the next OSPH file line
            if (scaleOnXi):
                if (file_line.startswith("  XI=")):
                    xi = float(file_line[5:20].replace('D', 'E'))
                    alam = 2.0 * math.pi * xi / vk
            else:
                if (file_line.startswith("  VK=")):
                    # we found VK, so we calculate lambda
                    # extract VK as a number from a string section, after
                    # replacing FORTRAN's D with E
                    vk = float(file_line[5:20].replace("D", "E"))
                    alam = 2.0 * math.pi / vk
            if ("SPHERE  1" in file_line):
                # we are in average section. We start a nested loop to
                # extract the average values
                found_averages = False
                while (not found_averages):
                    file_line = osph_file.readline()
                    if ("----- SCS ----- ABS ----- EXS ----- ALBEDS --" in file_line):
                        file_line = osph_file.readline()
                        # we now extract the values from string sections
                        scasm = float(file_line[1:15].replace("D", "E"))
                        abssm = float(file_line[17:30].replace("D", "E"))
                        extsm = float(file_line[32:45].replace("D", "E"))
                        # we can write the average values, similarly to fort.22
                        # note that \n puts a new line at the end of the string
                        output_line = "{0:.7E},{1:.7E},{2:.7E},{3:.7E}\n".format(alam, scasm, abssm, extsm)
                        if (out20 is not None): out20.write(output_line)
                        # we know that the asymmetry parameter and the radiation
                        # pressure forces will be after 5 more lines
                        for i in range(5):
                            file_line = osph_file.readline()
                        cosav = float(file_line[8:23].replace("D", "E"))
                        rapr = float(file_line[31:46].replace("D", "E"))
                        output_line = "{0:.7E},{1:.7E},{2:.7E}\n".format(alam, cosav, rapr)
                        if (out30 is not None): out30.write(output_line)
                        found_averages = True # terminate the inner loop
                # the averages were written. We look for the differential
                # section
                found_differentials = False
                while (not found_differentials):
                    for di in range(ndirs):
                        while ("JTH =" not in file_line):
                            file_line = osph_file.readline()
                        file_line = osph_file.readline()
                        tidg = float(file_line[7:17].replace("D", "E"))
                        pidg = float(file_line[24:34].replace("D", "E"))
                        tsdg = float(file_line[41:51].replace("D", "E"))
                        psdg = float(file_line[58:68].replace("D", "E"))
                        while ("SPHERE  1" not in file_line):
                            file_line = osph_file.readline()
                        if ("SPHERE  1" in file_line):
                            # we found SPHERE section. We know that the forces
                            # will be after 3 more lines
                            for i in range(3):
                                file_line = osph_file.readline()
                                # the following check is needed to parse C++ output
                                if ("INSERTION" in file_line):
                                    file_line = osph_file.readline()
                            fx = float(file_line[5:20].replace("D", "E"))
                            fy = float(file_line[25:40].replace("D", "E"))
                            fz = float(file_line[45:60].replace("D", "E"))
                            # we can write the differential values, similarly to fort.31
                            output_line = "{0:.7E},{1:.3E},{2:.3E},{3:.3E},{4:.3E},{5:.7E},{6:.7E},{7:.7E}\n".format(alam, tidg, tsdg, pidg, psdg, fx, fy, fz)
                            if (out40 is not None): out40.write(output_line)
                    found_differentials = True # terminate the inner loop
            # The parsing loop ends here

        if (out20 is not None): out20.close()
        if (out30 is not None): out30.close()
        if (out40 is not None): out40.close()
        osph_file.close() # close the OSPH file
    except Exception as ex:
        print(ex)
        errors += 1
    return errors

## \brief Print a command-line help summary.
def print_help():
    print("                                                  ")
    print("***               PARSE RESULTS                ***")
    print("                                                  ")
    print("Parse the results of a NPTMcode model calculation.")
    print("                                                                                          ")
    print("Usage: \"./parse_results.py --in INPUT --out OUTPUT [OPTIONS]\"                           ")
    print("                                                                                          ")
    print("Valid options are:                                                                        ")
    print("--in INPUT               File containing the results of the model calculation (mandatory).")
    print("--out OUTPUT             Root name for the output CSV data files (mandatory).             ")
    print("--app=[SPH|CLU|INCLU]    Application whose output needs to be parsed (mandatory).         ")
    print("--format=HDF5|LEGACY     Format of the result file to be parsed (autodetected by default).")
    print("--selection=[ALL|...]    Select the data to tabulate. Default is ALL. Optional filters can")
    print("                         be provided as a bracketed lists of the following values: ICS (i.")
    print("                         e. Integrated Cross Sections), DCS (i. e. Differential Cross Sec-")
    print("                         tions), IRP (i. e. Integrated Radiation Pressures), or DRP (i. e.")
    print("                         Differential Radiation Pressures)                                ")
    print("--help                   Print this help and exit.                                        ")
    print("                                                                                          ")

# ### PROGRAM EXECUTION ###
## \cond
res = main()
## \endcond
exit(res)
