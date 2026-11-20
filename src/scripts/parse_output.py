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
#  \ingroup python_tools
#  \brief Script to extract physical quantities from NPTMcode output files.
#
#  This script is intended to assist users in the inspection of results. It
#  extracts cross-sections, asymmetry parameters, radiation pressure forces
#  and torques, as functions of wavelength, from the output file of
#  `np_sphere`, `np_cluster` or `np_inclusion` and saves them in machine
#  readable CSV files.
#
#  The output file may be the legacy one (`c_OSPH`, `c_OCLU`, `c_OINCLU`,
#  or the output of the FORTRAN codes) or the HDF5 one (`c_OSPH.hd5`,
#  `c_OCLU.hd5`, `c_OINCLU.hd5`): both give the same CSV files. The format
#  is recognized from the file contents.
#
#  The CSV files are named after the root given with `--out`:
#  - `ROOT_ics1.csv`, `ROOT_ics2.csv`: integrated (ensemble average)
#    cross-sections in the polarization states -1 and +1;
#  - `ROOT_irp1.csv`, `ROOT_irp2.csv`: integrated asymmetry parameter and
#    radiation pressure cross-section;
#  - `ROOT_dcs1.csv`, `ROOT_dcs2.csv`: cross-sections for each incident and
#    scattered direction;
#  - `ROOT_drp1.csv`, `ROOT_drp2.csv`: asymmetry parameter, radiation pressure
#    forces and torques for each direction in which the code computes them
#    (every direction with ISAM < 0, the first scattered direction of each
#    incident direction otherwise).
#
#  For `np_sphere` the files are `ROOT_ics.csv`, `ROOT_irp.csv` (values of
#  the first sphere type) and `ROOT_drp.csv` (forces for each incident
#  direction).
#
#  The wavelength is 2 pi / VK, or 2 pi / (XI VK) when XI is a scale factor
#  for lengths (IDFC < 0 in the scatterer file).
#
#  The script execution requires python3; reading HDF5 files also requires
#  the `h5py` module.

import math
import re
import sys

from nptm_estimates import angle_count
from nptm_version import __version__

## \brief The programs whose output can be parsed, by application identifier.
PROGRAMS = {'SPH': 'np_sphere', 'CLU': 'np_cluster', 'INCLU': 'np_inclusion'}

## \brief The data selections accepted by `--selection`.
SELECTIONS = ('ALL', 'ICS', 'DCS', 'IRP', 'DRP')

## \brief Header of the integrated cross-section files.
ICS_HEADER = "Wavelength_m,ScaSec_m2,AbsSec_m2,ExtSec_m2"

## \brief Header of the integrated radiation pressure files.
IRP_HEADER = "Wavelength_m,CosAv,RaPr_m2"

## \brief Header of the differential cross-section files.
DCS_HEADER = "Wavelength_m,THi_deg,THs_deg,PHi_deg,PHs_deg,ScaSec_m2,AbsSec_m2,ExtSec_m2"

## \brief Header of the differential radiation pressure files of `np_cluster` and `np_inclusion`.
DRP_HEADER = (
    "Wavelength_m,THi_deg,THs_deg,PHi_deg,PHs_deg,CosAv,RaPr_m2,Fl_m2,Fr_m2,Fk_m2,"
    "Fx_m2,Fy_m2,Fz_m2,TQEl_m2,TQEr_m2,TQEk_m2,TQEx_m2,TQEy_m2,TQEz_m2,"
    "TQSl_m2,TQSr_m2,TQSk_m2,TQSx_m2,TQSy_m2,TQSz_m2"
)

## \brief Header of the differential radiation pressure file of `np_sphere`.
SPHERE_DRP_HEADER = "Wavelength_m,THi_deg,THs_deg,PHi_deg,PHs_deg,Fx_m2,Fy_m2,Fz_m2"

## \brief Labels of the legacy output for the columns of the differential
#  radiation pressure files of `np_cluster` and `np_inclusion`, after the
#  angles, in their order.
DRP_LABELS = (
    'COSAV', 'RAPRS', 'Fl', 'Fr', 'Fk', 'Fx', 'Fy', 'Fz',
    'TQEl', 'TQEr', 'TQEk', 'TQEx', 'TQEy', 'TQEz',
    'TQSl', 'TQSr', 'TQSk', 'TQSx', 'TQSy', 'TQSz'
)

## \brief Names of the HDF5 datasets that hold the columns of the
#  differential radiation pressure files, after the angles: a `%d` stands
#  for the polarization state (1 for -1, 2 for +1).
DRP_DATASETS = {
    'CLU': tuple("VEC_DIR_" + q + "C%d" for q in (
        'COSAV', 'RAPR', 'FL', 'FR', 'FK', 'FX', 'FY', 'FZ',
        'TQEL', 'TQER', 'TQEK', 'TQEX', 'TQEY', 'TQEZ',
        'TQSL', 'TQSR', 'TQSK', 'TQSX', 'TQSY', 'TQSZ')),
    'INCLU': tuple("VEC_DIR_" + q + "%d" for q in (
        'COSAV', 'RAPR', 'FL', 'FR', 'FK', 'FX', 'FY', 'FZ',
        'TQEL', 'TQER', 'TQEK', 'TQEX', 'TQEY', 'TQEZ',
        'TQSL', 'TQSR', 'TQSK', 'TQSX', 'TQSY', 'TQSZ'))
}

## \brief Names of the HDF5 datasets of the integrated values, by
#  application: scattering, absorption and extinction cross-sections,
#  asymmetry parameter and radiation pressure cross-section (`%d` as above).
AVERAGE_DATASETS = {
    'SPH': ('VEC_SCS', 'VEC_ABS', 'VEC_EXS', 'VEC_COSAV', 'VEC_RAPRS'),
    'CLU': ('VEC_SCC%d', 'VEC_ABC%d', 'VEC_EXC%d', 'VEC_COSAVC%d', 'VEC_RAPRC%d'),
    'INCLU': ('VEC_SCS%d', 'VEC_ABS%d', 'VEC_EXS%d', 'VEC_COSAV%d', 'VEC_RAPRS%d')
}

## \brief Names of the HDF5 datasets of the differential cross-sections (`%d` as above).
DCS_DATASETS = {
    'CLU': ('VEC_DIR_SCC%d', 'VEC_DIR_ABC%d', 'VEC_DIR_EXC%d'),
    'INCLU': ('VEC_DIR_SCS%d', 'VEC_DIR_ABS%d', 'VEC_DIR_EXS%d')
}

## \brief Names of the HDF5 datasets of the angle sweeps, in the order of
#  `direction_list()`.
ANGLE_DATASETS = (
    'TH_START', 'TH_STEP', 'TH_END', 'THS_START', 'THS_STEP', 'THS_END',
    'PH_START', 'PH_STEP', 'PH_END', 'PHS_START', 'PHS_STEP', 'PHS_END'
)

## \brief First bytes of an HDF5 file.
HDF5_SIGNATURE = b'\x89HDF\r\n\x1a\n'

## \brief A real number as written in the legacy files (also FORTRAN `D` exponents).
FLOAT_PATTERN = r'[-+]?(?:\d+\.\d*|\.\d+)(?:[DdEe][-+]?\d+)?'

## \brief A `LABEL=value` pair of the legacy files.
LABEL_VALUE = re.compile(r'([A-Za-z][A-Za-z0-9]*)=\s*(' + FLOAT_PATTERN + ')')

## \brief A line that starts a polarization state section (`LIN -1`, `CIRC  1`).
POLARIZATION_LINE = re.compile(r'^\s*(?:LIN|CIRC)\s+(-1|1)\s*$')

## \brief A line that starts the section of a sphere (`SPHERE  1`).
SPHERE_LINE = re.compile(r'^\s*SPHERE\s+(\d+)\s*$')

## \brief Error in the input data or in the arguments.
class ParseError(Exception):
    pass

## \brief Main execution code
#
# `main()` reads the arguments, parses the result file and writes the CSV
# files.
#
# \returns errors: `int` Number of errors (0 means successful run).
def main():
    try:
        config = parse_arguments(sys.argv[1:])
    except ParseError as ex:
        print("ERROR: %s" % ex)
        print("\nType \"parse_output.py --help\" to get more detailed help.")
        return 1
    if config['version_mode']:
        print("parse_output.py v%s." % __version__)
        return 0
    if config['help_mode']:
        print_help()
        return 0
    errors = 0
    if (config['result_file'] == ''):
        print("ERROR: no input file name given (missing \"--in FILE\" option).")
        errors += 1
    if (config['output_name'] == ''):
        print("ERROR: no output file name given (missing \"--out FILE\" option).")
        errors += 1
    if (config['application'] == ''):
        print("ERROR: no application identifier given (missing \"--app=APP\" option).")
        errors += 1
    if (errors > 0):
        return errors
    try:
        file_format = config['format']
        if (file_format == ''):
            file_format = detect_format(config['result_file'])
        if (file_format == 'HDF5'):
            tables = read_hdf5(config['result_file'], config['application'])
        else:
            tables = read_legacy(config['result_file'], config['application'])
        if (all(len(tables[name]) == 0 for name in tables)):
            raise ParseError(
                "no results found in %s: is it an output of %s?" % (
                    config['result_file'], PROGRAMS[config['application']]
                )
            )
        written = write_tables(tables, config['output_name'], config['selection'])
    except UnicodeDecodeError:
        print("ERROR: %s is not a legacy output file." % config['result_file'])
        return 1
    except (ParseError, OSError) as ex:
        print("ERROR: %s" % ex)
        return 1
    wavelengths = len(set(row[0] for name in tables for row in tables[name]))
    print(
        "%s: %d wavelengths read from %s (%s); written %s." % (
            "parse_output.py", wavelengths, config['result_file'], file_format,
            ", ".join(written) if len(written) > 0 else "no files"
        )
    )
    return 0

## \brief Parse the command line arguments.
#
#  The mandatory arguments are the name of the result file, the root name of
#  the CSV files and the application. The format is recognized from the file
#  contents, unless it is given.
#
#  \param args: `list(str)` The command line arguments, without the script name.
#  \returns config: `dict` A dictionary containing the script configuration.
def parse_arguments(args):
    config = {
        'result_file': '',
        'output_name': '',
        'application': '',
        'format': '',
        'selection': ['ALL'],
        'help_mode': False,
        'version_mode': False
    }
    options = {'--in': 'result_file', '--out': 'output_name', '--app': 'application',
               '--format': 'format', '--selection': 'selection'}
    i = 0
    while (i < len(args)):
        arg = args[i]
        name, sep, value = arg.partition('=')
        if (arg == "--help"):
            config['help_mode'] = True
        elif (arg == "--version"):
            config['version_mode'] = True
        elif (name in options):
            if (sep == ''):
                if (i + 1 >= len(args)):
                    raise ParseError("missing value of option \"%s\"" % name)
                i += 1
                value = args[i]
            config[options[name]] = value
        else:
            raise ParseError("unrecognized argument \"%s\"" % arg)
        i += 1
    config['application'] = config['application'].upper()
    if (config['application'] != '' and config['application'] not in PROGRAMS):
        raise ParseError(
            "unknown application \"%s\" (options are SPH|CLU|INCLU)" % config['application']
        )
    config['format'] = config['format'].upper()
    if (config['format'] not in ('', 'HDF5', 'LEGACY')):
        raise ParseError("unknown format \"%s\" (options are HDF5|LEGACY)" % config['format'])
    if (isinstance(config['selection'], str)):
        items = [s.strip().upper() for s in config['selection'].strip("[]").split(',')]
        for item in items:
            if (item not in SELECTIONS):
                raise ParseError(
                    "unknown selection \"%s\" (options are %s)" % (item, ", ".join(SELECTIONS))
                )
        config['selection'] = items
    return config

## \brief Recognize the format of a result file from its first bytes.
#
#  \param file_name: `str` The result file.
#  \returns format: `str` `HDF5` or `LEGACY`.
def detect_format(file_name):
    with open(file_name, "rb") as result_file:
        signature = result_file.read(len(HDF5_SIGNATURE))
    return 'HDF5' if (signature == HDF5_SIGNATURE) else 'LEGACY'

## \brief Empty tables for the CSV files of an application.
#
#  \param application: `str` The application identifier (`SPH`, `CLU` or `INCLU`).
#  \returns tables: `dict` Lists of rows, by file suffix (for instance `ics1`).
def new_tables(application):
    if (application == 'SPH'):
        return {'ics': [], 'irp': [], 'drp': []}
    return {
        'ics1': [], 'ics2': [], 'dcs1': [], 'dcs2': [],
        'irp1': [], 'irp2': [], 'drp1': [], 'drp2': []
    }

## \brief Name of the table of a quantity in a polarization state.
#
#  \param application: `str` The application identifier.
#  \param kind: `str` The quantity (`ics`, `irp`, `dcs` or `drp`).
#  \param state: `int` 1 for the state -1, 2 for the state +1 (ignored for `SPH`).
#  \returns name: `str` The table name (the suffix of the CSV file).
def table_name(application, kind, state):
    return kind if (application == 'SPH') else kind + str(state)

## \brief Wavelength of a scale.
#
#  \param vk: `float` Wave number (VK).
#  \param xi: `float` Scale value (XI).
#  \param xi_scales_lengths: `bool` Whether XI is a scale factor for lengths (IDFC < 0).
#  \returns wavelength: `float` The wavelength in meters.
def wavelength(vk, xi, xi_scales_lengths):
    if xi_scales_lengths:
        return 2.0 * math.pi / (xi * vk)
    return 2.0 * math.pi / vk

## \brief Scattered direction of the output files (`scattered_direction()` in the code).
#
#  With ISAM > 1 the scattered theta follows the incident one, with ISAM >= 1
#  the scattered phi is the incident one (plus 180 degrees across the poles).
#
#  \param isam: `int` The ISAM option.
#  \param th: `float` First incident theta.
#  \param ths: `float` First scattered theta.
#  \param tidg: `float` Incident theta.
#  \param pidg: `float` Incident phi.
#  \param ths_sweep: `float` Value of the scattered theta sweep.
#  \param phs_sweep: `float` Value of the scattered phi sweep.
#  \returns direction: `tuple(float, float)` Scattered theta and phi, in degrees.
def scattered_direction(isam, th, ths, tidg, pidg, ths_sweep, phs_sweep):
    tsdg = tidg + (ths - th) if (isam > 1) else ths_sweep
    phsph = 0.0
    if (isam >= 1):
        if (tsdg < 0.0 or tsdg > 180.0):
            phsph = 180.0
        if (tsdg < 0.0):
            tsdg = -tsdg
        if (tsdg > 180.0):
            tsdg = 360.0 - tsdg
    psdg = phs_sweep
    if (isam >= 1):
        psdg = pidg + phsph
        if (psdg >= 360.0):
            psdg -= 360.0
    return (tsdg, psdg)

## \brief Directions of the output, in the order of the code.
#
#  The order is incident theta, incident phi, scattered theta, scattered phi
#  (outermost to innermost). With ISAM > 1 there is one scattered theta, with
#  ISAM >= 1 one scattered phi.
#
#  \param angles: `list(float)` The 12 angles of the sweeps: incident theta
#  (start, step, end), scattered theta, incident phi, scattered phi.
#  \param isam: `int` The ISAM option.
#  \returns directions: `list(dict)` For each direction the angles `tidg`,
#  `tsdg`, `pidg`, `psdg` (degrees), the incident direction index `incident`
#  and whether it is the first scattered direction (`first`).
def direction_list(angles, isam):
    th, thstp, thlst, ths, thsstp, thslst, ph, phstp, phlst, phs, phsstp, phslst = angles
    num_theta = angle_count(th, thstp, thlst)
    num_thetas = 1 if (isam > 1) else angle_count(ths, thsstp, thslst)
    num_phi = angle_count(ph, phstp, phlst)
    num_phis = 1 if (isam >= 1) else angle_count(phs, phsstp, phslst)
    directions = []
    for jth in range(num_theta):
        for jph in range(num_phi):
            for jths in range(num_thetas):
                for jphs in range(num_phis):
                    tidg = th + jth * thstp
                    pidg = ph + jph * phstp
                    tsdg, psdg = scattered_direction(
                        isam, th, ths, tidg, pidg, ths + jths * thsstp, phs + jphs * phsstp
                    )
                    directions.append({
                        'tidg': tidg, 'tsdg': tsdg, 'pidg': pidg, 'psdg': psdg,
                        'incident': jth * num_phi + jph,
                        'first': (jths == 0 and jphs == 0)
                    })
    return directions

## \brief Read a dataset of an HDF5 result file as a one-dimensional array.
#
#  \param h5_file: `h5py.File` The open result file.
#  \param name: `str` The dataset name.
#  \param file_name: `str` The result file name (for the error message).
#  \param application: `str` The application identifier (for the error message).
#  \returns values: `numpy.ndarray` The values of the dataset.
def hdf5_dataset(h5_file, name, file_name, application):
    if (name not in h5_file):
        raise ParseError(
            "%s has no dataset %s: is it an output of %s?" % (
                file_name, name, PROGRAMS[application]
            )
        )
    return h5_file[name][()].reshape(-1)

## \brief Read the tables from an HDF5 result file.
#
#  \param file_name: `str` The result file.
#  \param application: `str` The application identifier.
#  \returns tables: `dict` The rows of the CSV files, by file suffix.
def read_hdf5(file_name, application):
    try:
        import h5py
    except ImportError:
        raise ParseError(
            "reading HDF5 files needs the Python module h5py (for instance: pip install h5py)."
        )
    tables = new_tables(application)
    with h5py.File(file_name, "r") as h5_file:
        data = lambda name: hdf5_dataset(h5_file, name, file_name, application)
        isam = int(data('ISAM')[0])
        idfc = int(data('IDFC')[0])
        nsph = int(data('NSPH')[0])
        if (application == 'SPH' and nsph != 1):
            raise ParseError(
                "the number of spheres is %d, only 1 is supported (is it an output of np_sphere?)" % nsph
            )
        configurations = int(data('NUM_CONF' if (application == 'SPH') else 'NCONF')[0])
        angles = [float(data(name)[0]) for name in ANGLE_DATASETS]
        directions = direction_list(angles, isam)
        ndirs = len(directions)
        num_incident = directions[-1]['incident'] + 1
        vec_vk = data('VEC_VK')
        vec_xi = data('VEC_XI')
        vec_ier = data('VEC_IER')
        vec_jxi = data('VEC_JXI')
        states = (0,) if (application == 'SPH') else (1, 2)
        averages = {}
        for state in states:
            averages[state] = [
                data(name % state if '%' in name else name)
                for name in AVERAGE_DATASETS[application]
            ]
        if (application == 'SPH'):
            forces = [data(name) for name in ('VEC_DIR_FX', 'VEC_DIR_FY', 'VEC_DIR_FZ')]
        else:
            differentials = {}
            pressures = {}
            for state in states:
                differentials[state] = [data(name % state) for name in DCS_DATASETS[application]]
                pressures[state] = [data(name % state) for name in DRP_DATASETS[application]]
        for jxi in range(len(vec_vk)):
            if (vec_ier[jxi] != 0):
                print(
                    "WARNING: scale %d stopped with error code %d: skipped." % (
                        int(vec_jxi[jxi]), int(vec_ier[jxi])
                    )
                )
                continue
            alam = wavelength(float(vec_vk[jxi]), float(vec_xi[jxi]), idfc < 0)
            for state in states:
                # In np_sphere the first configuration is the first sphere type.
                ai = jxi * configurations if (application == 'SPH') else jxi
                values = [float(v[ai]) for v in averages[state]]
                tables[table_name(application, 'ics', state)].append([alam] + values[0:3])
                tables[table_name(application, 'irp', state)].append([alam] + values[3:5])
            for di in range(ndirs):
                d = directions[di]
                angle_values = [d['tidg'], d['tsdg'], d['pidg'], d['psdg']]
                if (application == 'SPH'):
                    if d['first']:
                        fi = (jxi * num_incident + d['incident']) * nsph
                        tables['drp'].append(
                            [alam] + angle_values + [float(v[fi]) for v in forces]
                        )
                    continue
                k = jxi * ndirs + di
                for state in states:
                    tables['dcs%d' % state].append(
                        [alam] + angle_values + [float(v[k]) for v in differentials[state]]
                    )
                    if (isam < 0 or d['first']):
                        tables['drp%d' % state].append(
                            [alam] + angle_values + [float(v[k]) for v in pressures[state]]
                        )
    return tables

## \brief Convert a number of a legacy file.
#
#  \param text: `str` The number, possibly with a FORTRAN `D` exponent.
#  \returns value: `float` The number.
def to_float(text):
    return float(text.replace('D', 'E').replace('d', 'e'))

## \brief Read the tables from a legacy result file.
#
#  The values are found by their labels and by the section headers, as the
#  code and the FORTRAN programs write them: the ensemble averages
#  (`CLUSTER (ENSEMBLE AVERAGE, ...)`, `ENSEMBLE AVERAGE, ...` or the first
#  sphere of `np_sphere`), and for each direction (`JTH = ...`) the values of
#  the particle (`CLUSTER`, `SINGLE SCATTERER` or the first sphere), in the
#  polarization states -1 and +1. Lines starting with `INSERTION` are
#  skipped.
#
#  \param file_name: `str` The result file.
#  \param application: `str` The application identifier.
#  \returns tables: `dict` The rows of the CSV files, by file suffix.
def read_legacy(file_name, application):
    tables = new_tables(application)
    scale_vk = None      # VK, when XI is a scale factor for lengths
    alam = None          # wavelength of the current scale
    in_directions = False
    target = None        # 'avg', 'dir' or None (values not extracted)
    state = 0            # polarization state: 1 for -1, 2 for +1
    direction = None     # angles of the current direction
    record = None        # values of the current radiation pressure row
    expect_cross_sections = False
    with open(file_name, "r") as result_file:
        for line_number, line in enumerate(result_file, start=1):
            stripped = line.strip()
            if stripped.startswith("INSERTION"):
                continue
            where = "%s, line %d" % (file_name, line_number)
            if expect_cross_sections:
                expect_cross_sections = False
                numbers = re.findall(FLOAT_PATTERN, line)
                if (len(numbers) < 3):
                    raise ParseError("%s: cross-sections expected." % where)
                values = [to_float(n) for n in numbers[0:3]]
                if (alam is None):
                    raise ParseError("%s: no wavelength before the values." % where)
                if (target == 'avg'):
                    tables[table_name(application, 'ics', state)].append([alam] + values)
                elif (target == 'dir'):
                    tables['dcs%d' % state].append([alam] + direction + values)
                continue
            if (application == 'SPH' and line_number == 2):
                nsph = int(stripped.split()[0])
                if (nsph != 1):
                    raise ParseError(
                        "the number of spheres is %d, only 1 is supported (is it an output of np_sphere?)" % nsph
                    )
            labels = dict(LABEL_VALUE.findall(line))
            if ("SCALE FACTOR FOR LENGTHS" in line):
                scale_vk = to_float(labels['VK'])
            elif ("========== JXI =" in line):
                alam = None
                in_directions = False
                target = None
            elif ((stripped.startswith("VK=") or stripped.startswith("XI=")) and not in_directions):
                if (scale_vk is not None):
                    alam = wavelength(scale_vk, to_float(labels['XI']), True)
                else:
                    alam = wavelength(to_float(labels['VK']), 1.0, False)
            elif ("JTH =" in line):
                in_directions = True
                target = None
                direction = None
            elif ('TIDG' in labels):
                direction = [
                    to_float(labels['TIDG']), to_float(labels['TSDG']),
                    to_float(labels['PIDG']), to_float(labels['PSDG'])
                ]
            elif (POLARIZATION_LINE.match(line) and target is not None):
                state = 1 if (POLARIZATION_LINE.match(line).group(1) == '-1') else 2
            elif stripped.startswith("----- SC"):
                if (target is not None):
                    if (application != 'SPH' and state == 0):
                        raise ParseError("%s: values without a polarization state." % where)
                    if (target == 'dir' and direction is None):
                        raise ParseError("%s: values without a direction." % where)
                    expect_cross_sections = True
            elif (SPHERE_LINE.match(line)):
                first = (SPHERE_LINE.match(line).group(1) == '1')
                if (application == 'SPH' and first):
                    target = 'dir' if in_directions else 'avg'
                else:
                    target = None
            elif (application == 'CLU' and stripped.startswith("CLUSTER (ENSEMBLE AVERAGE")):
                target, state = 'avg', 0
            elif (application == 'CLU' and stripped == "CLUSTER"):
                target, state = 'dir', 0
            elif (application == 'INCLU' and stripped.startswith("ENSEMBLE AVERAGE")):
                target, state = 'avg', 0
            elif (application == 'INCLU' and stripped == "SINGLE SCATTERER"):
                target, state = 'dir', 0
            elif (target == 'avg' and 'COSAV' in labels):
                if (alam is None):
                    raise ParseError("%s: no wavelength before the values." % where)
                tables[table_name(application, 'irp', state)].append(
                    [alam, to_float(labels['COSAV']), to_float(labels['RAPRS'])]
                )
            elif (target == 'dir' and application == 'SPH' and 'Fx' in labels):
                if (direction is None):
                    raise ParseError("%s: forces without a direction." % where)
                tables['drp'].append(
                    [alam] + direction + [to_float(labels[q]) for q in ('Fx', 'Fy', 'Fz')]
                )
            elif (target == 'dir' and application != 'SPH'):
                if ('COSAV' in labels):
                    record = {}
                if (record is not None):
                    record.update(labels)
                    if ('TQSz' in labels):
                        missing = [q for q in DRP_LABELS if q not in record]
                        if (len(missing) > 0):
                            raise ParseError("%s: missing %s." % (where, ", ".join(missing)))
                        tables['drp%d' % state].append(
                            [alam] + direction + [to_float(record[q]) for q in DRP_LABELS]
                        )
                        record = None
    return tables

## \brief Write the selected tables as CSV files.
#
#  \param tables: `dict` The rows of the CSV files, by file suffix.
#  \param root_name: `str` The root name of the CSV files.
#  \param selection: `list(str)` The selected quantities (`ALL`, `ICS`, `DCS`, `IRP`, `DRP`).
#  \returns written: `list(str)` The names of the written files.
def write_tables(tables, root_name, selection):
    written = []
    for name in tables:
        kind = name[0:3]
        if ('ALL' not in selection and kind.upper() not in selection):
            continue
        if (kind == 'ics'):
            header = ICS_HEADER
        elif (kind == 'irp'):
            header = IRP_HEADER
        elif (kind == 'dcs'):
            header = DCS_HEADER
        else:
            header = SPHERE_DRP_HEADER if (name == 'drp') else DRP_HEADER
        file_name = root_name + "_" + name + ".csv"
        with open(file_name, "w") as csv_file:
            csv_file.write(header + "\n")
            for row in tables[name]:
                if (kind in ('dcs', 'drp')):
                    fields = (
                        ["{0:.7E}".format(row[0])]
                        + ["{0:.3E}".format(v) for v in row[1:5]]
                        + ["{0:.7E}".format(v) for v in row[5:]]
                    )
                else:
                    fields = ["{0:.7E}".format(v) for v in row]
                csv_file.write(",".join(fields) + "\n")
        written.append(file_name)
    return written

## \brief Print a command-line help summary.
def print_help():
    print("                                                  ")
    print("***               PARSE RESULTS                ***")
    print("                                                  ")
    print("Parse the results of a NPTMcode model calculation.")
    print("                                                                                          ")
    print("Usage: \"./parse_output.py --in INPUT --out OUTPUT --app=APP [OPTIONS]\"                  ")
    print("                                                                                          ")
    print("Valid options are:                                                                        ")
    print("--in INPUT               File containing the results of the model calculation, legacy or  ")
    print("                         HDF5 (mandatory).                                                ")
    print("--out OUTPUT             Root name for the output CSV data files (mandatory).             ")
    print("--app=[SPH|CLU|INCLU]    Application whose output needs to be parsed (mandatory).         ")
    print("--format=HDF5|LEGACY     Format of the result file to be parsed (autodetected by default).")
    print("--selection=[ALL|...]    Select the data to tabulate. Default is ALL. Optional filters can")
    print("                         be provided as a bracketed lists of the following values: ICS (i.")
    print("                         e. Integrated Cross Sections), DCS (i. e. Differential Cross Sec-")
    print("                         tions), IRP (i. e. Integrated Radiation Pressures), or DRP (i. e.")
    print("                         Differential Radiation Pressures)                                ")
    print("--help                   Print this help and exit.                                        ")
    print("--version                Print the version and exit.                                      ")
    print("                                                                                          ")
    print("Reading HDF5 files requires the Python module h5py.                                       ")
    print("                                                                                          ")

## \cond
if __name__ == "__main__":
    sys.exit(main())
## \endcond
