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

## @package nptm_estimates
#  \ingroup python_tools
#  \brief Particle radius and memory estimates, as computed by the C++ code.
#
#  This module reproduces, in python, the quantities that the programs log
#  before a calculation, so that the scripts can report them without running
#  the code:
#
#  - the radius of the smallest sphere centred in the origin of the
#    coordinates that contains the particle
#    (`ScattererConfiguration::get_particle_radius()`), which sets the
#    external expansion order with dynamic orders;
#  - the host memory of the iteration data, of the inversion workspace and
#    of the output (the `get_size()` functions of the C++ classes, and
#    `get_inversion_workspace()` in `utils.cpp`);
#  - the GPU memory of an inversion (`get_device_inversion_size()`) and the
#    margins of the GPU memory check (`resources.h`).
#
#  The C++ functions are the reference: a change in one of them must be
#  copied here. The sizes of the C++ class instances are those of an x86_64
#  Linux build.
#
#  The module is imported by `guess_resources.py` and `model_maker.py`.

import math

## \cond
GIB = 1024.0 * 1024.0 * 1024.0
DCOMPLEX = 16
DOUBLE = 8
INT = 4
LONG = 8
POINTER = 8
SHORT = 2
BOOL = 1
# sizeof() of the C++ classes whose instance is counted (x86_64 Linux).
SIZEOF_SPHERE_OUTPUT_INFO = 600
SIZEOF_CLUSTER_OUTPUT_INFO = 1536
SIZEOF_INCLUSION_OUTPUT_INFO = 1184
SIZEOF_PARTICLE_DESCRIPTOR_SPHERE = 976
SIZEOF_SPHERE_ITERATION_DATA = 440
## \endcond

## \brief GPU memory set aside for each MPI process that uses a GPU, in bytes.
#
#  Same as `NPTM_GPU_PROCESS_OVERHEAD` in `resources.h`.
GPU_PROCESS_OVERHEAD = 256 * 1024 * 1024

## \brief Fraction of the free GPU memory that a calculation may use.
#
#  Same as `NPTM_GPU_USABLE_FRACTION` in `resources.h`.
GPU_USABLE_FRACTION = 0.9

## \brief Columns of the LU workspaces counted on the GPU (largest MAGMA block size).
LU_BLOCK_COLUMNS = 512

## \brief Inversion back ends that a build can use.
BACKENDS = ("magma", "cublas", "lapack", "lucin")

## \brief Inversion modes accepted by the `INVERSION` option.
INVERSION_MODES = ("LU", "GESV", "RBT", "SVD")

## \brief Radius of the smallest sphere centred in the origin that contains the particle.
#
#  This is the radius used by the code (`get_particle_radius()`) to choose
#  the external expansion order: the external field is expanded around the
#  origin of the coordinates, so the radius is max(|r_i| + a_i) over the
#  spheres. For an inclusion, it is the radius of the external sphere.
#
#  \param radii: `list(float)` Radius of each sphere, in meters.
#  \param vec_x: `list(float)` x coordinate of each sphere, in meters.
#  \param vec_y: `list(float)` y coordinate of each sphere, in meters.
#  \param vec_z: `list(float)` z coordinate of each sphere, in meters.
#  \param external_radius: `float` Radius of the external sphere of an inclusion (optional; None for other particles).
#  \return radius: `float` The encircling radius, in meters.
def particle_radius(radii, vec_x, vec_y, vec_z, external_radius=None):
    if (external_radius is not None):
        return external_radius
    result = 0.0
    for i in range(len(radii)):
        dist = math.sqrt(vec_x[i] * vec_x[i] + vec_y[i] * vec_y[i] + vec_z[i] * vec_z[i]) + radii[i]
        if (dist > result):
            result = dist
    return result

## \brief Number of values of an angle sweep (`angle_sweep_count()` in the code).
#
#  \param start: `float` First value.
#  \param step: `float` Step (0 for a single value).
#  \param end: `float` Last value.
#  \return count: `int` Number of values.
def angle_count(start, step, end):
    if (step == 0.0):
        return 1
    return 1 + int((end - start) / step + 1.0e-3)

## \brief Numbers of incident and scattering directions.
#
#  With ISAM > 1 there is one scattered theta and one scattered phi, with
#  ISAM = 1 one scattered phi.
#
#  \param angles: `list(float)` The 12 angles of the geometry file: incident
#  theta (start, step, end), scattered theta, incident phi, scattered phi.
#  \param isam: `int` The ISAM option of the geometry file (optional, default 0).
#  \return counts: `tuple(int, int, int, int, int)` Number of incident
#  theta, scattered theta, incident phi and scattered phi values, and number
#  of directions (their product).
def direction_counts(angles, isam=0):
    num_theta = angle_count(angles[0], angles[1], angles[2])
    num_thetas = 1 if isam > 1 else angle_count(angles[3], angles[4], angles[5])
    num_phi = angle_count(angles[6], angles[7], angles[8])
    num_phis = 1 if isam >= 1 else angle_count(angles[9], angles[10], angles[11])
    return (num_theta, num_thetas, num_phi, num_phis, num_theta * num_thetas * num_phi * num_phis)

## \brief Layers of the sphere types.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return layers: `tuple(int, int)` Total number of layers and largest number of layers of a type.
def _layers(model):
    num_layers = 0
    max_layers = 1
    for nli in range(model['configurations']):
        nl = model['nshl'][nli]
        if (nli == 0 and model['app'] == "INCLUSION"):
            nl += 1
        num_layers += nl
        if (nl > max_layers):
            max_layers = nl
    return (num_layers, max_layers)

## \brief Order-dependent sizes shared by the cluster and inclusion descriptors.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return sizes: `dict` The derived sizes, named as in the C++ code.
def _orders(model):
    nsph = model['nsph']
    li = model['li']
    le = model['le']
    s = {}
    s['lm'] = li if li > le else le
    s['nlim'] = li * (li + 2)
    s['nlem'] = le * (le + 2)
    s['nlemt'] = 2 * s['nlem']
    s['ncou'] = nsph * (nsph - 1)
    s['litpo'] = li + li + 1
    s['litpos'] = s['litpo'] * s['litpo']
    s['lmtpo'] = li + le + 1
    s['lmtpos'] = s['lmtpo'] * s['lmtpo']
    lm = s['lm']
    s['nv3j'] = (lm * (lm + 1) * (2 * lm + 7)) // 6
    s['ndi'] = nsph * s['nlim']
    s['ndit'] = 2 * nsph * s['nlim']
    s['nllt'] = 2 * nsph * li * (li + 2) if s['nlemt'] == 0 else s['nlemt']
    s['ndm'] = 2 * (nsph * s['nlim'] + s['nlem'])
    s['nhspo'] = 2 * max(model['npnt'], model['npntts']) - 1
    return s

## \brief Host memory of the common part of the cluster and inclusion descriptors.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \param s: `dict` Derived sizes (from `_orders()`).
#  \return bytes: `int` Size in bytes.
def _descriptor_common(model, s):
    nsph = model['nsph']
    li = model['li']
    nconf = model['configurations']
    num_layers, max_layers = _layers(model)
    result = SHORT * 2 + INT * 16 + LONG * 20 + DOUBLE
    result += DCOMPLEX * (2 * li * nsph + 16 + 2 * s['nhspo'] + nsph + max_layers + 1)
    result += DCOMPLEX * (2 * li + nconf)
    result += DOUBLE * (num_layers + nconf + 4 * nsph)
    result += INT * (nsph + nconf)
    return result

## \brief Host memory of the "cluster and inclusion" section of the descriptors.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \param s: `dict` Derived sizes (from `_orders()`).
#  \param pointers: `int` Number of pointers of the section.
#  \return bytes: `int` Size in bytes.
def _descriptor_cluster_inclusion(model, s, pointers):
    nsph = model['nsph']
    lm = s['lm']
    result = INT * 29 + LONG * pointers + DCOMPLEX
    result += DCOMPLEX * (
        4 * s['nllt'] + s['nlemt'] * s['nlemt'] + 36 + s['ncou'] * s['litpo']
        + nsph * s['lmtpo'] + s['ncou'] * s['litpos'] + nsph * s['lmtpos']
    )
    result += INT * ((lm + 1) * lm)
    result += DOUBLE * (8 + s['nv3j'] + s['lmtpo'])
    return result

## \brief Host memory of `ParticleDescriptorCluster`.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return bytes: `int` Size in bytes.
def cluster_descriptor_size(model):
    s = _orders(model)
    nsph = model['nsph']
    result = _descriptor_common(model, s)
    result += LONG * 12 + DCOMPLEX * (22 * nsph) + DOUBLE * (7 * nsph) + LONG * (4 * nsph)
    result += LONG * 9 + DCOMPLEX + DOUBLE * 3
    result += DCOMPLEX * (20 + 2 * s['ndi'] * s['nlem'] + s['ndit'] * s['nlemt'])
    result += LONG * (2 + 2 * s['ndit'])
    result += _descriptor_cluster_inclusion(model, s, 31)
    return result

## \brief Host memory of `ParticleDescriptorInclusion`.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return bytes: `int` Size in bytes.
def inclusion_descriptor_size(model):
    s = _orders(model)
    result = _descriptor_common(model, s)
    result += _descriptor_cluster_inclusion(model, s, 25)
    result += DCOMPLEX * (8 * model['le'] + s['nlemt'] * s['ndm'])
    result += LONG * s['nlemt']
    return result

## \brief Host memory of `ParticleDescriptorSphere` (base class arrays included).
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return bytes: `int` Size in bytes.
def sphere_descriptor_size(model):
    nsph = model['nsph']
    li = model['lm'] if nsph == 1 else model['li']
    nconf = model['configurations']
    num_layers, max_layers = _layers(model)
    nhspo = 2 * max(model['npnt'], model['npntts']) - 1
    nllt = 2 * nsph * li * (li + 2)
    # ParticleDescriptor::get_size()
    base = 2 * DCOMPLEX * li * nsph + 2 * POINTER * li + DOUBLE * num_layers
    base += POINTER * nconf + DCOMPLEX * 16 + 4 * DOUBLE * nsph + INT * nsph
    base += DOUBLE * nconf + INT * nconf + 2 * DCOMPLEX * nhspo + DCOMPLEX * nsph
    base += DCOMPLEX * (max_layers + 1)
    result = SIZEOF_PARTICLE_DESCRIPTOR_SPHERE + base
    result += DCOMPLEX * 4 * nsph + DCOMPLEX * 16 * nsph + DCOMPLEX * nllt * 4
    result += POINTER * nllt + DCOMPLEX * nsph + 7 * DOUBLE * nsph
    result += POINTER * nsph + POINTER * nsph + POINTER * 2 * nsph
    return result

## \brief Host memory of the iteration data of one thread.
#
#  `ClusterIterationData`, `InclusionIterationData` or `SphereIterationData`,
#  with their particle descriptor. For clusters and inclusions they include
#  the matrix to invert.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return bytes: `int` Size in bytes.
def iteration_data_size(model):
    nsph = model['nsph']
    if (model['app'] == "SPHERE"):
        lm = model['lm']
        result = SIZEOF_SPHERE_ITERATION_DATA + sphere_descriptor_size(model)
        result += 2 * DOUBLE + DOUBLE * nsph + 11 * DOUBLE * 3 + 2 * DOUBLE * 16
        result += 2 * POINTER * 4 + 2 * DCOMPLEX * 2 * nsph + 2 * DOUBLE * 2 * nsph
        result += 2 * POINTER * 2 + 2 * POINTER * 2
        result += DOUBLE * lm * 12 + POINTER * lm + POINTER * 3 * lm + POINTER * 2 * 3 * lm
        return result
    s = _orders(model)
    lm = s['lm']
    if (model['app'] == "CLUSTER"):
        n = s['ndit']
        result = LONG * 37 + DOUBLE * 9 + DCOMPLEX + INT * 5 + BOOL
        descriptor = cluster_descriptor_size(model)
    else:
        n = s['ndm']
        result = LONG * 37 + DOUBLE * 10 + DCOMPLEX + INT * 6 + BOOL
        descriptor = inclusion_descriptor_size(model)
    result += LONG * (44 + 6 * lm + n)
    result += DOUBLE * (132 + 5 * nsph + 12 * lm)
    result += DCOMPLEX * (24 + 4 * nsph + n * n)
    return result + descriptor

## \brief Host memory of the output data.
#
#  `ClusterOutputInfo`, `InclusionOutputInfo` or `SphereOutputInfo` for
#  `scales` scales.
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \param scales: `int` Number of scales (optional; 0 for all the scales of the model).
#  \return bytes: `int` Size in bytes.
def output_size(model, scales=0):
    nsph = model['nsph']
    nconf = model['configurations']
    xib = model['nxi'] if scales == 0 else scales
    nt, nts, nph, nphs, ndirs = direction_counts(model['angles'], model.get('isam', 0))
    angles = DOUBLE * (nt + nph + nts + nphs) + 5 * DOUBLE * ndirs + 2 * DOUBLE * 3 * ndirs
    if (model['app'] == "SPHERE"):
        result = SIZEOF_SPHERE_OUTPUT_INFO + INT * xib + SHORT * xib + 2 * DOUBLE * xib
        result += 17 * DOUBLE * nconf * xib + 2 * DCOMPLEX * nconf * xib + angles
        result += 4 * DCOMPLEX * nsph * ndirs * xib + 3 * DOUBLE * nsph * nt * nph * xib
        result += 2 * DOUBLE * 16 * nsph * ndirs * xib
        if (nsph != 1):
            result += DCOMPLEX * xib + 3 * DOUBLE * xib
    elif (model['app'] == "CLUSTER"):
        result = SIZEOF_CLUSTER_OUTPUT_INFO + 3 * DOUBLE * nsph + INT * xib + SHORT * xib
        result += 37 * DOUBLE * xib + 17 * DOUBLE * xib * nconf + 2 * DCOMPLEX * xib * nconf
        result += 5 * DCOMPLEX * xib + angles
        result += 4 * DCOMPLEX * ndirs * nconf * xib + 2 * DOUBLE * 16 * ndirs * nconf * xib
        result += 12 * DCOMPLEX * ndirs * xib + 66 * DOUBLE * ndirs * xib
        result += 2 * DOUBLE * 16 * ndirs * xib
    else:
        result = SIZEOF_INCLUSION_OUTPUT_INFO + 3 * DOUBLE * nsph + INT * xib + SHORT * xib
        result += 28 * DOUBLE * xib + DOUBLE * xib * (nconf + 1) + DCOMPLEX * xib * (nconf + 1)
        result += 4 * DCOMPLEX * xib + angles
        result += 60 * DOUBLE * ndirs * xib + 8 * DCOMPLEX * ndirs * xib
        result += 2 * DOUBLE * 16 * ndirs * xib
    return result

## \brief Order of the matrix to invert (0 for a single sphere).
#
#  \param model: `dict` Model description (see `memory_needs()`).
#  \return order: `int` Number of rows (and columns) of the matrix.
def matrix_order(model):
    if (model['app'] == "SPHERE"):
        return 0
    s = _orders(model)
    return s['ndit'] if model['app'] == "CLUSTER" else s['ndm']

## \brief Host workspace of a CPU inversion (LAPACK, or `lucin()` without LAPACK).
#
#  \param n: `int` Order of the matrix.
#  \param mode: `str` Inversion mode (`LU`, `GESV`, `RBT` or `SVD`).
#  \param refine: `bool` Whether the iterative refinement is enabled.
#  \param lapack: `bool` Whether the build uses LAPACK.
#  \param int_bytes: `int` Size of the LAPACK integers (8 with ILP64, 4 otherwise).
#  \return bytes: `int` Size in bytes.
def cpu_inversion_workspace(n, mode, refine, lapack, int_bytes):
    matrix = DCOMPLEX * n * n
    if (not lapack):
        return 2 * DOUBLE * n
    result = int_bytes * n
    if (mode == "GESV"):
        result += matrix
    if (refine):
        result += 4 * matrix + DCOMPLEX * n
    return result

## \brief Host workspace of one inversion, besides the matrix (`get_inversion_workspace()`).
#
#  \param n: `int` Order of the matrix.
#  \param mode: `str` Inversion mode (`LU`, `GESV`, `RBT` or `SVD`).
#  \param refine: `bool` Whether the iterative refinement is enabled.
#  \param backend: `str` Inversion back end of the build (one of `BACKENDS`).
#  \param lapack: `bool` Whether the build uses LAPACK (for the CPU inversion).
#  \param int_bytes: `int` Size of the LAPACK and MAGMA integers.
#  \param on_gpu: `bool` Whether the GPU inverts the matrix (MAGMA and cuBLAS builds).
#  \return bytes: `int` Size in bytes.
def inversion_workspace(n, mode, refine, backend, lapack, int_bytes, on_gpu=True):
    matrix = DCOMPLEX * n * n
    if (backend == "magma" and on_gpu):
        result = int_bytes * n
        if (mode == "GESV"):
            result += matrix
        elif (mode == "RBT"):
            result += 2 * matrix
        elif (mode == "SVD"):
            result += 3 * matrix + DOUBLE * (5 * n * n + 7 * n)
        if (refine):
            result += DCOMPLEX * n
        return result
    if (backend == "cublas" and on_gpu):
        return DCOMPLEX
    # The CPU inversion of a MAGMA or cuBLAS build uses LAPACK if available.
    use_lapack = (backend == "lapack") or (backend in ("magma", "cublas") and lapack)
    return cpu_inversion_workspace(n, mode, refine, use_lapack, int_bytes)

## \brief GPU memory of one inversion, matrix included (`get_device_inversion_size()`).
#
#  \param n: `int` Order of the matrix.
#  \param mode: `str` Inversion mode (`LU`, `GESV`, `RBT` or `SVD`).
#  \param refine: `bool` Whether the iterative refinement is enabled.
#  \param backend: `str` Inversion back end of the build (one of `BACKENDS`).
#  \return bytes: `int` Size in bytes (0 for CPU back ends).
def device_inversion_size(n, mode, refine, backend):
    matrix = DCOMPLEX * n * n
    vector = DCOMPLEX * n
    if (backend == "magma"):
        if (mode == "GESV"):
            return 5 * matrix + vector if refine else 2 * matrix + LU_BLOCK_COLUMNS * vector
        if (mode == "RBT"):
            return (5 if refine else 2) * matrix + 4 * vector
        if (mode == "SVD"):
            return 6 * matrix + vector if refine else 3 * matrix
        result = matrix + 2 * LU_BLOCK_COLUMNS * vector
        if (refine and 4 * matrix + vector > result):
            result = 4 * matrix + vector
        return result
    if (backend == "cublas"):
        result = 2 * matrix + LU_BLOCK_COLUMNS * vector + 8 * n
        if (refine):
            result += 2 * matrix
        return result
    return 0

## \brief Memory needs of a calculation, as the programs compute them.
#
#  The model description is a dictionary with the keys:
#  - `app`: `"SPHERE"`, `"CLUSTER"` or `"INCLUSION"`;
#  - `nsph`, `lm` (single sphere), `li`, `le`, `npnt`, `npntts`, `nxi`;
#  - `configurations` and `nshl` (layers of each sphere type);
#  - `angles`: the 12 angles of the geometry file, and `isam` (default 0);
#  - `refine` (`bool`), `inv_mode` (`str`), `device_policy` (`str`).
#
#  \param model: `dict` Model description.
#  \param backend: `str` Inversion back end of the build (one of `BACKENDS`, optional, default `lapack`).
#  \param threads: `int` OpenMP threads of each MPI process (optional, default 1).
#  \param lapack: `bool` Whether a MAGMA or cuBLAS build also uses LAPACK (optional, default True).
#  \param int_bytes: `int` Size of the LAPACK and MAGMA integers (optional, default 8 for ILP64).
#  \return needs: `dict` Sizes in bytes: `iteration_data`, `workspace`,
#  `output`, `scale_output` (per thread or per process, as named),
#  `process_0` and `other_process` (host memory of an MPI process),
#  `matrix_order`, `matrix`, `device_per_thread` and `device_per_process`
#  (GPU memory, with `GPU_PROCESS_OVERHEAD`; 0 without GPU).
def memory_needs(model, backend="lapack", threads=1, lapack=True, int_bytes=8):
    needs = {}
    n = matrix_order(model)
    mode = model.get('inv_mode', "LU")
    refine = bool(model.get('refine', False))
    needs['matrix_order'] = n
    needs['matrix'] = DCOMPLEX * n * n
    needs['iteration_data'] = iteration_data_size(model)
    needs['output'] = output_size(model)
    needs['scale_output'] = output_size(model, 1)
    if (n > 0):
        workspace = inversion_workspace(n, mode, refine, backend, lapack, int_bytes)
        if (model.get('device_policy', "STOP") == "USE_CPU_FALLBACK" and backend in ("magma", "cublas")):
            cpu = inversion_workspace(n, mode, refine, backend, lapack, int_bytes, False)
            if (cpu > workspace):
                workspace = cpu
        device = device_inversion_size(n, mode, refine, backend)
    else:
        workspace = 0
        device = 0
    needs['workspace'] = workspace
    per_thread = needs['iteration_data'] + workspace
    needs['process_0'] = needs['output'] + threads * per_thread + (threads - 1) * needs['scale_output']
    needs['other_process'] = threads * per_thread + (2 * threads - 1) * needs['scale_output']
    needs['device_per_thread'] = device
    needs['device_per_process'] = threads * device + GPU_PROCESS_OVERHEAD if device > 0 else 0
    return needs
