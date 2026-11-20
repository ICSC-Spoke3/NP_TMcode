#!/bin/bash

# Default configuration settings
ARFLAGS="-rs"
BUILDFORTRAN="auto"
CUBLAS="auto"
CXX_DBG=" -ggdb"
CXX_OPT=3
DEBUGFLAGS=""
ENABLE_ILP64="yes"
FC_OPT=3
FC_DBG=" -ggdb"
FFLAGS=""
INCLUDEFLAGS=""
LAPACK="auto"
LDFLAGS=""
LIBMODE="static"
MAGMA="auto"
MAGMA_INVERT_FLAGS=""
NVTXFLAGS=""
OMPMODE="auto"
OFFLOAD="auto"
OFFLOADFLAGS=""
OFFLOADLDFLAGS=""
OMPLDFLAGS=""
# Link only the libraries that the programs use (the default of g++ on some
# systems, not of clang): an unused library can otherwise provide symbols to
# other libraries, e.g. 64-bit-integer LAPACK routines to a 32-bit MAGMA.
ASNEEDED="-Wl,--as-needed "
OFFLOAD_ARCH="auto"
OFFLOAD_CUDA_PATH="auto"
OFFLOAD_TARGET="auto"
MAGMA_BACKEND="none"
# End of default configuration settings

# Function declarations
function guess_cxx {
    # Guess the name of the C++ compiler
    result=$(which mpicxx 2>/dev/null)
    if [ "x$result" = "x" ]; then
	result=$(which g++ 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result=$(which clang 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result=$(which icpx 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result="none"
    fi
    echo $result
}

function guess_fc {
    # Guess the name of the FORTRAN compiler
    result=$(which mpif90 2>/dev/null)
    if [ "x$result" = "x" ]; then
	result=$(which gfortran 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result=$(which f77 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result=$(which flang 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result=$(which flang-new 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result=$(which ifx 2>/dev/null)
    fi
    if [ "x$result" = "x" ]; then
	result="none"
    fi
    echo $result
}

function print_help {
    echo "NPTMcode configuration script.                                                   "
    echo "                                                                                 "
    echo "Use this script to detect the proper build configuration for your system.        "
    echo "Valid options are:                                                               "
    echo "                                                                                 "
    echo "--enable-fortran                  Enable FORTRAN compilation (DEFAULT).          "
    echo "--disable-fortran                 Disable FORTRAN compilation.                   "
    echo "--enable-gdb                      Enable GNU debugger (DEFAULT).                 "
    echo "--disable-gdb                     Disable GNU debugger.                          "
    echo "--enable-ilp64                    Enable 64-bit integers (DEFAULT).              "
    echo "--disable-ilp64                   Disable 64-bit integers.                       "
    echo "--enable-offload                  Enable GPU offloading (DEFAULT).               "
    echo "--disable-offload                 Disable GPU offloading.                        "
    echo "--with-offload-target=TARGET      Offload to NVIDIA (nvptx) or AMD (amdgcn) GPUs."
    echo "                                  (AMD support is experimental: not yet tested.) "
    echo "                                  Default: the GPUs of MAGMA (CUDA or HIP), else "
    echo "                                  the GPUs found by nvidia-smi or amdgpu-arch.   "
    echo "--with-offload-arch=ARCH          Compile offloaded code for GPU architecture    "
    echo "                                  ARCH (e.g. sm_89, gfx90a). Default: detected   "
    echo "                                  (nvidia-smi, amdgpu-arch or rocminfo), else    "
    echo "                                  sm_70 / gfx90a with g++, native with clang.    "
    echo "--with-offload-cuda=PATH          CUDA installation used by clang for offloading "
    echo "                                  (default: detected).                           "
    echo "--enable-openmp                   Enable OpenMP multi-threading (DEFAULT).       "
    echo "--disable-openmp                  Disable OpenMP multi-threading.                "
    echo "--enable-debug=FEATURE            Enable debug output of specified feature.      "
    echo "--enable-nvtx                     Enable NVTX profiling tools.                   "
    echo "--disable-nvtx                    Disable NVTX profiling tools (DEFAULT).        "
    echo "--enable-optimize=OPT             Use OPT-level compiler optimization.           "
    echo "--enable-shared                   Use shared libraries (default is static).      "
    echo "--disable-shared                  Use static libraries (DEFAULT).                "
    echo "--help                            Print this help and exit.                      "
    echo "--with-cublas                     Use cuBLAS and cuSOLVER (DEFAULT).             "
    echo "--without-cublas                  Disable cuBLAS.                                "
    echo "--with-fflags=FFLAGS              Use specified FORTRAN compiler flags instead of"
    echo "                                  detected ones.                                 "
    echo "--with-hdf5=HDF5_PATH             Use specified HDF5 distribution.               "
    echo "--with-include=[PATH]             Use additional include paths (can be given more"
    echo "                                  than once).                                    "
    echo "--with-lapack                     Use LAPACK (DEFAULT).                          "
    echo "--without-lapack                  Disable LAPACK.                                "
    echo "--with-ldflags=[PATH]             Use additional linker LDFLAGS (can be given    "
    echo "                                  more than once).                               "
    echo "--with-magma=[MAGMA]              Use specified MAGMA distribution (DEFAULT).    "
    echo "--without-magma                   Disable MAGMA.                                 "
    echo "--enable-magma-invert=OPTION      Use specified MAGMA matrix inversion function  "
    echo "                                  (options = auto|getri|gesv|rbt, default getri)."
    echo "                                                                                 "
    echo "Some influential environment variables are:                                      "
    echo "                                                                                 "
    echo "AR                                Static library archive packager.               "
    echo "CUDA_HOME (or CUDA_PATH)          Root folder of the CUDA installation to use    "
    echo "                                  (default: the newest one found by pkg-config)."
    echo "CUDAFLAGS                         CUDA compilation flags.                        "
    echo "CUDALDFLAGS                       CUDA libraries linking flags.                  "
    echo "CXX                               C++ compiler.                                  "
    echo "CXXFLAGS                          C++ compilation flags.                         "
    echo "CXXLDFLAGS                        C++ linker flags.                              "
    echo "FC                                FORTRAN compiler.                              "
    echo "HDF5_INCLUDE                      Path to the HDF5 header files.                 "
    echo "HDF5_LIB                          Path to the HDF5 library files.                "
    echo "MAGMA_HOME                        Root folder of MAGMA distribution.             "
    echo "ROCM_PATH                         Root folder of the ROCm installation, for a    "
    echo "                                  MAGMA built on HIP and for AMD offloading      "
    echo "                                  (default: /opt/rocm). Experimental.            "
    echo "OMPI_CXX, MPICH_CXX               C++ compiler used by the Open MPI or MPICH     "
    echo "                                  wrapper (e.g. clang++ for clang offloading).   "
    echo "                                  Recorded in Makefile.in for the build.         "
}

function test_legacy_fortran {
    cat > conf_test_fortran.f <<EOF
      PROGRAM CONF_TEST_FORTRAN
      I=2
      J=I-I
      END
EOF
    $1 -std=legacy conf_test_fortran.f -o conf_test_fortran > /dev/null 2>>error.log
    result=$?
    rm conf_test_fortran.f
    if [ "x$result" = "x0" ]; then
	rm conf_test_fortran
    fi
    echo $result
}
# End of function declarations

# Script execution
num_args=${#@}
declare -a args=("$@")

# Argument parsing section
for arg in "${args[@]}"
do
    cut_arg=$(echo $arg | cut -d '=' -f1)
    if [ "x$cut_arg" = "x--help" ]; then
	print_help | less
	exit 0
    elif [ "x$cut_arg" = "x--enable-gdb" ]; then
	FC_DBG=" -ggdb"
	CXX_DBG=" -ggdb"
    elif [ "x$cut_arg" = "x--disable-gdb" ]; then
	FC_DBG=""
	CXX_DBG=""
    elif [ "x$cut_arg" = "x--enable-ilp64" ]; then
	ENABLE_ILP64="yes"
    elif [ "x$cut_arg" = "x--disable-ilp64" ]; then
	ENABLE_ILP64="no"
    elif [ "x$cut_arg" = "x--enable-fortran" ]; then
	BUILDFORTRAN="yes"
    elif [ "x$cut_arg" = "x--disable-fortran" ]; then
	FC="none"
	BUILDFORTRAN="no"
    elif [ "x$cut_arg" = "x--enable-offload" ]; then
	OFFLOAD="auto"
    elif [ "x$cut_arg" = "x--disable-offload" ]; then
	OFFLOAD="no"
    elif [ "x$cut_arg" = "x--with-offload-target" ]; then
	OFFLOAD_TARGET=$(echo $arg | cut -d '=' -f2)
	if [ "x$OFFLOAD_TARGET" != "xnvptx" -a "x$OFFLOAD_TARGET" != "xamdgcn" ]; then
	    echo "ERROR: unknown offload target \"$OFFLOAD_TARGET\" (use nvptx or amdgcn)."
	    exit 1
	fi
    elif [ "x$cut_arg" = "x--with-offload-arch" ]; then
	OFFLOAD_ARCH=$(echo $arg | cut -d '=' -f2)
    elif [ "x$cut_arg" = "x--with-offload-cuda" ]; then
	OFFLOAD_CUDA_PATH=$(echo $arg | cut -d '=' -f2)
    elif [ "x$cut_arg" = "x--enable-openmp" ]; then
	OMPMODE="yes"
    elif [ "x$cut_arg" = "x--disable-openmp" ]; then
	OMPMODE="no"
    elif [ "x$cut_arg" = "x--enable-debug" ]; then
	dbg_feature=$(echo $arg | cut -d '=' -f2)
	if [ "x$dbg_feature" = "x" ]; then
	    echo "ERROR: no debug feature specified!"
	    echo "ERROR: no debug feature specified!" >>configure.log
	    exit 1
	else
	    if [ "x$DEBUGFLAGS" = "x" ]; then
		DEBUGFLAGS=" -DDEBUG_${dbg_feature}"
	    else
		DEBUGFLAGS="${DEBUGFLAGS} -DDEBUG_${dbg_feature}"
	    fi
	fi
    elif [ "x$cut_arg" = "x--enable-nvtx" ]; then
	NVTXFLAGS=" -DUSE_NVTX"
    elif [ "x$cut_arg" = "x--disable-nvtx" ]; then
	NVTXFLAGS=""
    elif [ "x$cut_arg" = "x--enable-optimize" ]; then
	opt_level=$(echo $arg | cut -d '=' -f2)
	if [ $opt_level -lt 0 -o $opt_level -gt 3 ]; then
	    echo "ERROR: invalid optimization level $opt_level"
	    echo "ERROR: invalid optimization level $opt_level" >>configure.log
	    exit 1
	fi
	FC_OPT=$opt_level
	CXX_OPT=$opt_level
    elif [ "x$cut_arg" = "x--enable-shared" ]; then
	LIBMODE="shared"
    elif [ "x$cut_arg" = "x--disable-shared" ]; then
	LIBMODE="static"
    elif [ "x$cut_arg" = "x--with-cublas" ]; then
	CUBLAS="yes"
    elif [ "x$cut_arg" = "x--without-cublas" ]; then
	CUBLAS="no"
    elif [ "x$cut_arg" = "x--with-fflags" ]; then
	custom_flags=$(echo $arg | cut -d '=' -f2)
	if [ "x$custom_flags" != "x" ]; then
	    num_characters=${#custom_flags}
	    if [ $num_characters -gt 4 ]; then
		if [ "x${custom_flags:$num_characters-4:4}" = "x-std" ]; then
		    std_flags=$(echo $arg | cut -d '=' -f3)
		    custom_flags="${custom_flags}=${std_flags}"
		fi
	    fi
	    FFLAGS=$custom_flags
	fi
    elif [ "x$cut_arg" = "x--with-hdf5" ]; then
	HDF5_HOME=$(echo $arg | cut -d '=' -f2)
	if [ "x${HDF5_HOME}" = "x" ]; then
	    echo "ERROR: option --with-hdf5 requires a valid HDF5 path."
	    echo "ERROR: option --with-hdf5 requires a valid HDF5 path." >>configure.log
	    exit 1
	fi
    elif [ "x$cut_arg" = "x--with-include" ]; then
	user_flags=$(echo $arg | cut -d '=' -f2)
	INCLUDEFLAGS="$INCLUDEFLAGS $user_flags"
    elif [ "x$cut_arg" = "x--with-lapack" ]; then
	LAPACK="yes"
    elif [ "x$cut_arg" = "x--without-lapack" ]; then
	LAPACK="no"
    elif [ "x$cut_arg" = "x--with-ldflags" ]; then
	user_flags=$(echo $arg | cut -d '=' -f2)
	LDFLAGS="$LDFLAGS $user_flags"
    elif [ "x$cut_arg" = "x--with-magma" ]; then
	MAGMA="yes"
	MAGMA_HOME=$(echo $arg | cut -d '=' -f2)
    elif [ "x$cut_arg" = "x--without-magma" ]; then
	MAGMA="no"
    elif [ "x$cut_arg" = "x--enable-magma-invert" ]; then
	MAGMA_INVERT_CHOICE=$(echo $arg | cut -d '=' -f2)
	if [ "x$MAGMA_INVERT_CHOICE" = "xgetri" -o "x$MAGMA_INVERT_CHOICE" = "xauto" ]; then
	    MAGMA_INVERT_FLAGS=""
	elif [ "x$MAGMA_INVERT_CHOICE" = "xgesv" ]; then
	    MAGMA_INVERT_FLAGS=" -DUSE_ZGESV_GPU"
	elif [ "x$MAGMA_INVERT_CHOICE" = "xrbt" ]; then
	    MAGMA_INVERT_FLAGS=" -DUSE_ZGESV_RBT"
	else
	    echo "ERROR: unrecognized --enable-magma-invert option \"$MAGMA_INVERT_CHOICE\""
	    echo "ERROR: unrecognized --enable-magma-invert option \"$MAGMA_INVERT_CHOICE\"" >>configure.log
	    exit 1
	fi
    else
	echo "ERROR: unrecognized argument \"$arg\""
	echo "ERROR: unrecognized argument \"$arg\"" >>configure.log
	exit 1
    fi
done
# End of argument parsing section

# Configuration logic
echo "NPtm_code configuration"
echo "NPtm_code configuration" > configure.log
# The version of the package is the "version" key of CITATION.cff, in the
# root folder of the repository; the programs get it as NPTM_VERSION_ID.
echo -n "configure: checking the code version... "
echo -n "configure: checking the code version... " >>configure.log
CITATION_FILE="$(cd "$(dirname "$0")/.." && pwd)/CITATION.cff"
NPTM_VERSION=""
if [ -f "$CITATION_FILE" ]; then
    NPTM_VERSION=$(sed -n 's/^version:[[:space:]]*["'"'"']\{0,1\}\([^"'"'"'[:space:]]*\).*/\1/p' "$CITATION_FILE" | head -n 1)
fi
echo "$NPTM_VERSION" | grep -q '^[A-Za-z0-9][A-Za-z0-9._+-]*$'
if [ "x$?" = "x0" ]; then
    echo "$NPTM_VERSION"
    echo "$NPTM_VERSION" >>configure.log
    VERSIONFLAGS=" -DNPTM_VERSION_ID=${NPTM_VERSION}"
else
    echo "unknown"
    echo "unknown" >>configure.log
    echo "WARNING: no valid \"version:\" key in $CITATION_FILE; the programs will report an unknown version."
    echo "WARNING: no valid \"version:\" key in $CITATION_FILE; the programs will report an unknown version." >>configure.log
    VERSIONFLAGS=""
fi
echo -n "" > error.log
# Check for AR
echo -n "configure: checking for ar... "
echo -n "configure: checking for ar... " >>configure.log
if [ "x$AR" = "x" ]; then
    AR=$(which ar)
    if [ "x$AR" = "x" ]; then
	AR="ar"
    fi
fi
which $AR >/dev/null 2>>error.log
result=$?
if [ "x$result" != "x0" ]; then
    echo "none"
    echo "none" >>configure.log
    if [ "x$LIBMODE" = "xstatic" ]; then
	echo "ERROR: $AR not found!"
	echo "ERROR: $AR not found!" >>configure.log
	exit 2
    fi
else
    echo $AR
    echo $AR >>configure.log
fi
# Check FORTRAN compiler (if required)
if [ "x$FC" = "x" ]; then
    echo -n "configure: checking for FORTRAN compiler... "
    echo -n "configure: checking for FORTRAN compiler... " >>configure.log
    FC=$(guess_fc)
    echo $FC
    echo $FC >>configure.log
fi
if [ "x$FC" = "xnone" ]; then
    if [ "x$BUILDFORTRAN" = "xyes" ]; then
	echo "ERROR: FORTRAN compilation was requested, but no FORTRAN compiler found."
	echo "ERROR: FORTRAN compilation was requested, but no FORTRAN compiler found." >>configure.log
	exit 2
    else
	BUILDFORTRAN=""
    fi
else
    if [ "x$BUILDFORTRAN" = "xauto" ]; then
	BUILDFORTRAN="yes"
    fi
    echo -n "configure: checking whether $FC supports -ggdb... "
    echo -n "configure: checking whether $FC supports -ggdb... " >>configure.log
    cat > test_fortran.f <<EOF
      PROGRAM CONF_TEST_FORTRAN
      I=2
      J=I-I
      END
EOF
    $FC -ggdb test_fortran.f -o test_fortran > /dev/null 2>>error.log
    result=$?
    if [ "x$result" = "x0" ]; then
	echo "yes"
	echo "yes" >>configure.log
	rm test_fortran.f test_fortran
    else
	echo "no"
	echo "no" >>configure.log
	rm test_fortran.f
	FC_DBG=""
    fi
    echo -n "configure: checking whether $FC supports legacy... "
    echo -n "configure: checking whether $FC supports legacy... " >>configure.log
    result=$(test_legacy_fortran $FC)
    if [ "x$result" = "x0" ]; then
	FCFLAGS="-O${FC_OPT}${FC_DBG} -std=legacy"
	echo "yes"
	echo "yes" >>configure.log
    else
	FCFLAGS="-O$FC_OPT${FC_DBG}"
	echo "no"
	echo "WARNING: FORTRAN compiler does not support legacy flag."
	echo "no" >>configure.log
	echo "WARNING: FORTRAN compiler does not support legacy flag." >>configure.log
    fi
fi # End of FORTRAN compiler check

# Check C++ compiler (mandatory)
if [ "x$CXX" = "x" ]; then
    echo -n "configure: checking for C++ compiler... "
    echo -n "configure: checking for C++ compiler... " >>configure.log
    CXX=$(guess_cxx)
    echo $CXX
    echo $CXX >>configure.log
else
    echo "configure: using $CXX as C++ compiler."
    echo "configure: using $CXX as C++ compiler." >>configure.log
fi
if [ "x$CXX" = "xnone" ]; then
    echo "ERROR: no C++ compiler found!"
    echo "ERROR: no C++ compiler found!" >>configure.log
    exit 2
fi
CLANGFLAGS=""
$CXX --version | grep "clang" > /dev/null
result=$?
if [ "x$result" = "x0" ]; then
    CLANGFLAGS=" -stdlib=libstdc++"
fi
echo -n "configure: checking whether $CXX works... "
echo -n "configure: checking whether $CXX works... " >>configure.log
cat > test_compiler.cpp <<EOF
int main() {
  int i = -1;
  int j = 1;
  return (i + j);
}
EOF
$CXX $CLANGFLAGS test_compiler.cpp -o test_compiler > /dev/null 2>>error.log
result=$?
if [ "x$result" = "x0" ]; then
    echo "yes"
    echo "yes" >>configure.log
else
    echo "no"
    echo "ERROR: $CXX is not a working C++ compiler!"
    echo "no" >>configure.log
    echo "ERROR: $CXX is not a working C++ compiler!" >>configure.log
    exit 2
fi
echo -n "configure: checking whether $CXX supports -ggdb... "
echo -n "configure: checking whether $CXX supports -ggdb... " >>configure.log
$CXX $CLANGFLAGS -ggdb test_compiler.cpp -o test_compiler > /dev/null 2>>error.log
result=$?
if [ "x$result" = "x0" ]; then
    echo "yes"
    echo "yes" >>configure.log
    rm test_compiler.cpp test_compiler
else
    echo "no"
    echo "no" >>configure.log
    rm test_compiler.cpp
    CXX_DBG=""
fi
echo -n "configure: checking whether $CXX is a GNU compiler... "
echo -n "configure: checking whether $CXX is a GNU compiler... " >>configure.log
$CXX --version | grep "g++" > /dev/null
result=$?
if [ "x$result" = "x0" ]; then
    echo "yes"
    echo "yes" >>configure.log
else
    echo "no"
    echo "no" >>configure.log
fi
echo -n "configure: checking whether $CXX is a MPI compiler... "
echo -n "configure: checking whether $CXX is a MPI compiler... " >>configure.log
cat > test_compiler.cpp <<EOF
# include <mpi.h>
int main() {
  int result;
#ifdef MPI_VERSION
  result = 0;
#else
  result = 1;
#endif
  return result;
}
EOF
$CXX test_compiler.cpp -o test_compiler > /dev/null 2>>error.log
result=$?
if [ "x$result" = "x0" ]; then
    ./test_compiler
    result=$?
    if [ "x$result" = "x0" ]; then
	echo "yes"
	echo "yes" >>configure.log
	MPIFLAGS=" -DUSE_MPI"
    else
	echo "no"
	echo "no" >>configure.log
	MPIFLAGS=""
    fi
    rm test_compiler test_compiler.cpp
else
    echo "no"
    echo "no" >>configure.log
    MPIFLAGS=""
    rm test_compiler.cpp
fi
if [ "x$OMPMODE" != "xno" ]; then
    echo -n "configure: checking whether $CXX supports OpenMP... "
    echo -n "configure: checking whether $CXX supports OpenMP... " >>configure.log
    cat > test_compiler.cpp <<EOF
#include <omp.h>
int main() {
  int result = 0;
  if (omp_get_num_threads() < 1) result = 1;
  return result;
}
EOF
    $CXX -fopenmp test_compiler.cpp -o test_compiler > /dev/null 2>>error.log
    result=$?
    if [ "x$result" = "x0" ]; then
	./test_compiler
	result=$?
	if [ "x$result" = "x0" ]; then
	    echo "yes"
	    echo "yes" >>configure.log
	    OMPFLAGS=" -fopenmp -DUSE_OMP"
	else
	    echo "no"
	    echo "no" >>configure.log
	    OMP_FLAGS=""
	    OFFLOAD="no"
	    if [ "x$OMPMODE" = "xyes" ]; then
		echo "ERROR: OpenMP was requested, but it is not supported."
		echo "ERROR: OpenMP was requested, but it is not supported." >>configure.log
		rm test_compiler test_compiler.cpp
		exit 2
	    fi
	fi
	rm test_compiler test_compiler.cpp
    else
	echo "no"
	echo "no" >>configure.log
	OMP_FLAGS=""
	OFFLOAD="no"
	if [ "x$OMPMODE" = "xyes" ]; then
	    echo "ERROR: OpenMP was requested, but it is not supported."
	    echo "ERROR: OpenMP was requested, but it is not supported." >>configure.log
	    rm test_compiler.cpp
	    exit 2
	fi
	rm test_compiler.cpp
    fi
fi
if [ "x$OMPFLAGS" != "x" ]; then
    # clang implements OpenMP reductions of complex variables with 16-byte
    # atomic operations, which are library calls provided by libatomic.
    echo -n "configure: checking whether OpenMP complex reductions need libatomic... "
    echo -n "configure: checking whether OpenMP complex reductions need libatomic... " >>configure.log
    cat > test_compiler.cpp <<EOF
int main() {
  __complex__ double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < 100; i++) sum += 1.0;
  return (__real__ sum == 100.0) ? 0 : 1;
}
EOF
    $CXX${CLANGFLAGS} -fopenmp test_compiler.cpp -o test_compiler > /dev/null 2>>error.log
    if [ "x$?" = "x0" ]; then
	echo "no"
	echo "no" >>configure.log
    else
	$CXX${CLANGFLAGS} -fopenmp test_compiler.cpp -o test_compiler -latomic > /dev/null 2>>error.log
	if [ "x$?" = "x0" ]; then
	    echo "yes"
	    echo "yes" >>configure.log
	    OMPLDFLAGS=" -latomic"
	else
	    echo "unknown"
	    echo "unknown" >>configure.log
	fi
    fi
    rm -f test_compiler test_compiler.cpp
fi
# End of C++ compiler check
# Check HDF5
echo -n "configure: checking for HDF5 header flags..."
echo -n "configure: checking for HDF5 header flags..." >>configure.log
if [ "x$HDF5_HOME" != "x" ]; then
    HDF5_INCLUDE="$HDF5_HOME/include"
    HDF5_LIB="$HDF5_HOME/lib"
fi
if [ "x${HDF5_INCLUDE}${HDF5_LIB}" = "x" ]; then
    pkg-config --version > /dev/null
    use_pkg_config=$?
    if [ "x$use_pkg_config" = "x0" ]; then
	declare -a pkg_array=$(pkg-config --list-all | grep hdf5-serial)
	for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep hdf5-serial > /dev/null
	result=$?
	if [ "x$result" = "x0" ]; then
	    cflags=$(pkg-config --cflags-only-I hdf5-serial)
	    HDF5FLAGS=$(echo " ${cflags}")
	    ldflags=$(pkg-config --libs hdf5-serial)
	    HDF5LDFLAGS=$(echo "${ldflags}")
	else
	    declare -a pkg_array=$(pkg-config --list-all | grep hdf5)
	    for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep hdf5 > /dev/null
	    result=$?
	    if [ "x$result" = "x0" ]; then
		cflags=$(pkg-config --cflags-only-I hdf5)
		HDF5FLAGS=$(echo " ${cflags}")
		ldflags=$(pkg-config --libs hdf5)
		HDF5LDFLAGS=$(echo "${ldflags}")
	    fi
	fi
    else
	export -p | grep HDF5_ROOT > /dev/null
	result=$?
	if [ "x$result" = "x0" ]; then
	    if [ "x$HDF5FLAGS" = "x" ]; then
		HDF5FLAGS=" -I${HDF5_ROOT}/include"
	    fi
	    if [ "x$HDF5LDFLAGS" = "x" ]; then
		HDF5LDFLAGS="-L${HDF5_ROOT}/lib -lhdf5"
	    fi
	fi
	export -p | grep HDF5_DIR > /dev/null
	result=$?
	if [ "x$result" = "x0" ]; then
	    if [ "x$HDF5FLAGS" = "x" ]; then
		HDF5FLAGS=" -I${HDF5_DIR}/include"
	    fi
	    if [ "x$HDF5LDFLAGS" = "x" ]; then
		HDF5LDFLAGS="-L${HDF5_DIR}/lib -lhdf5"
	    fi
	fi
	if [ "x$HDF5FLAGS" = "x" ]; then
	    HDF5FLAGS=" -I/usr/include/hdf5/serial"
	fi
	if [ "x$HDF5LDFLAGS" = "x" ]; then
	    HDF5LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5"
	fi
    fi
else
    HDF5FLAGS=" -I${HDF5_INCLUDE}"
    HDF5LDFLAGS="-L${HDF5_LIB} -lhdf5"
fi
if [ "x$HDF5FLAGS" = "x" ]; then
    echo "not found."
    echo "ERROR: HDF5 headers not found!"
    echo "not found." >>configure.log
    echo "ERROR: HDF5 headers not found!" >>configure.log
    exit 2
else
    echo "$HDF5FLAGS"
    echo "$HDF5FLAGS" >>configure.log
fi
# End of HDF5 check
# LAPACK checks
if [ "x$LAPACK" != "xno" ]; then
    echo -n "configure: checking for LAPACK... "
    echo -n "configure: checking for LAPACK... " >>configure.log
    if [ "x$ENABLE_ILP64" = "xyes" ]; then
	# 64-bit indices are enabled
	LAPACK_ILP64_FLAG="-DLAPACK_ILP64 -DUSE_ILP64"
	LAPACK_ILP64_LDSPEC="_ilp64"
	LAPACK_LDSPEC="64"
	MKL_BUILD="mkl-dynamic-ilp64-gomp"
    else
	# 64-bit indices are disabled
	LAPACK_ILP64_FLAG=""
	LAPACK_ILP64_LDSPEC="_lp64"
	LAPACK_LDSPEC=""
	MKL_BUILD="mkl-dynamic-lp64-gomp"
    fi # end of 64-bit decision tree
    BLASFLAGS=""
    BLASLDFLAGS=""
    pkg-config --version > /dev/null
    use_pkg_config=$?
    if [ "x$use_pkg_config" = "x0" ]; then
	# pkg-config is available
	declare -a pkg_array=$(pkg-config --list-all | grep ${MKL_BUILD})
	for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep ${MKL_BUILD} > /dev/null
	result=$?
	if [ "x$result" = "x0" ]; then
            # MKL was found
            MKL_INCLUDE=$(pkg-config --cflags-only-I ${MKL_BUILD})
            LAPACKFLAGS=" -DUSE_LAPACK -DUSE_MKL ${LAPACK_ILP64_FLAG} ${MKL_INCLUDE}"
            LAPACKLDFLAGS="$(pkg-config --libs ${MKL_BUILD})"
	    if [ "x$CLANGFLAGS" != "x" ]; then
		# clang uses its own OpenMP runtime (libomp), not libgomp: link
		# the Intel threading layer of MKL, which uses the same entry
		# points as libomp, and not libgomp. Two OpenMP runtimes in one
		# process make omp_get_thread_num() return 0 in every thread.
		# Without OpenMP, nothing provides those entry points: use the
		# sequential layer.
		MKL_THREAD_LIB="-lmkl_intel_thread"
		if [ "x$OMPFLAGS" = "x" ]; then
		    MKL_THREAD_LIB="-lmkl_sequential"
		fi
		LAPACKLDFLAGS=$(echo " $LAPACKLDFLAGS " | sed -e "s/ -lmkl_gnu_thread / ${MKL_THREAD_LIB} /" -e 's/ -lgomp / /')
	    fi
	    echo "MKL"
	    echo "MKL" >>configure.log
      else
          # MKL was not found, so configuration searches for BLAS
          declare -a pkg_array=$(pkg-config --list-all | grep blas${LAPACK_LDSPEC})
          for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep blas${LAPACK_LDSPEC} > /dev/null
          result=$?
          if [ "x$result" = "x0" ]; then
              # BLAS was found
	      BLASLDFLAGS=$(pkg-config --libs blas${LAPACK_LDSPEC})
	  else
              declare -a pkg_array=$(pkg-config --list-all | grep openblas${LAPACK_LDSPEC})
              for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep openblas${LAPACK_LDSPEC} > /dev/null
              result=$?
              if [ "x$result" = "x0" ]; then
		  # OPENBLAS was found
		  BLASLDFLAGS=$(pkg-config --libs openblas${LAPACK_LDSPEC})
	      fi
          fi # end of BLAS decision tree
          # search for LAPACKe
          declare -a pkg_array=$(pkg-config --list-all | grep lapacke${LAPACK_LDSPEC})
          for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep lapacke${LAPACK_LDSPEC} > /dev/null
          result=$?
          if [ "x$result" = "x0" ]; then
              # LAPACKe was found
              LAPACK_INCLUDE=$(pkg-config --cflags-only-I lapacke${LAPACK_LDSPEC})
              LAPACKFLAGS=" -DUSE_LAPACK ${LAPACK_ILP64_FLAG} ${LAPACK_INCLUDE} ${BLASFLAGS}"
              LAPACKLDFLAGS="$(pkg-config --libs lapacke${LAPACK_LDSPEC}) ${BLASLDFLAGS}"
	      echo "lapacke"
	      echo "lapacke" >>configure.log
          fi # end of LAPACKe decision tree
          if [ "x${LAPACKFLAGS}${LAPACKLDFLAGS}" = "x" ]; then
              # LAPACKe was not found, so configuration searches for LAPACK
              declare -a pkg_array=$(pkg-config --list-all | grep lapack${LAPACK_LDSPEC})
              for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep lapack${LAPACK_LDSPEC} > /dev/null
              result=$?
              if [ "x$result" = "x0" ]; then
		  # LAPACK was found
		  LAPACK_INCLUDE=$(pkg-config --cflags-only-I lapack${LAPACK_LDSPEC})
		  LAPACKFLAGS=" -DUSE_LAPACK ${LAPACK_ILP64_FLAG} ${LAPACK_INCLUDE} ${BLASFLAGS}"
		  LAPACKLDFLAGS="$(pkg-config --libs lapack${LAPACK_LDSPEC}) ${BLASLDFLAGS}"
		  echo "LAPACK"
		  echo "LAPACK" >>configure.log
              fi # end of LAPACK decision tree
          fi # end of LAPACKe decision tree
	fi # end of MKL decision tree
    else
	# pkg-config is not available
	export -p | grep MKL > /dev/null
	MKL_DEF=$?
	if [ "x$MKL_DEF" = "x0" ]; then
            LAPACKFLAGS=" -DUSE_LAPACK -DUSE_MKL ${LAPACK_ILP64_FLAG} -I${MKLROOT}/include"
            if [ "x$CLANGFLAGS" != "x" ]; then
		# clang: Intel threading layer with clang's libomp (see above).
		MKL_THREAD_LIB="-lmkl_intel_thread"
		if [ "x$OMPFLAGS" = "x" ]; then
		    MKL_THREAD_LIB="-lmkl_sequential"
		fi
		LAPACKLDFLAGS=" -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel${LAPACK_ILP64_LDSPEC} ${MKL_THREAD_LIB} -lmkl_core -lpthread -lm -ldl"
	    else
		LAPACKLDFLAGS=" -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel${LAPACK_ILP64_LDSPEC} -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
	    fi
	    echo "MKL"
	    echo "MKL" >>configure.log
	else
            if [ -f /usr/include/lapacke.h ]; then
		LAPACKFLAGS=" -DUSE_LAPACK ${LAPACK_ILP64_FLAG} ${BLASFLAGS}"
		LAPACKLDFLAGS=" -llapacke${LAPACK_LDSPEC} ${BLASLDFLAGS}"
		echo "lapacke"
		echo "lapacke" >>configure.log
            fi
	fi
    fi
    if [ "x$LAPACKFLAGS" = "x" ]; then
	echo "no"
	echo "no" >>configure.log
	if [ "x$LAPACK" = "xyes" ]; then
	    echo "ERROR: LAPACK was required, but no LAPACK was found."
	    echo "ERROR: LAPACK was required, but no LAPACK was found." >>configure.log
	fi
    fi
fi
# End of LAPACK checks
# cuBLAS checks
if [ "x$CUBLAS" != "xno" ]; then
    echo -n "configure: checking for cuBLAS... "
    echo -n "configure: checking for cuBLAS... " >>configure.log
    pkg-config --version > /dev/null
    use_pkg_config=$?
    # A CUDA installation named by CUDA_HOME (or CUDA_PATH) takes precedence
    # over the ones that pkg-config lists, so that all the CUDA flags come
    # from a single version when several are installed.
    cuda_home="${CUDA_HOME:-${CUDA_PATH}}"
    if [ "x${CUDAFLAGS}${CUDALDFLAGS}" = "x" -a "x$cuda_home" != "x" ]; then
	if [ -f "${cuda_home}/include/cuda_runtime.h" ]; then
	    cuda_libdir="${cuda_home}/lib64"
	    if [ ! -d "$cuda_libdir" ]; then
		cuda_libdir="${cuda_home}/lib"
	    fi
	    CUDAFLAGS="-I${cuda_home}/include"
	    CUDALDFLAGS="-L${cuda_libdir} -lcublas -lcudart"
	fi
    fi
    if [ "x${CUDAFLAGS}${CUDALDFLAGS}" = "x" ]; then
	if [ "x$use_pkg_config" = "x0" ]; then
            # pkg-config is available
            declare -a pkg_array=$(pkg-config --list-all | grep cublas)
            for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cublas > /dev/null
            result=$?
            if [ "x$result" = "x0" ]; then
		# CUBLAS detected
		cuda_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cublas | sort -V | tail -n 1)
    		CUDAFLAGS=$(pkg-config --cflags ${cuda_pkg})
    		CUDALDFLAGS=$(pkg-config --libs ${cuda_pkg})
	    else
		# CUBLAS not detected
		CUBLAS="no"
		CUDAFLAGS=""
		CUDALDFLAGS=""
            fi # end of CUBLAS runtime decision tree
            declare -a pkg_array=$(pkg-config --list-all | grep cudart)
            for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart > /dev/null
            result=$?
            if [ "x$result" = "x0" ]; then
		# CUDA runtime detected
		cuda_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart | sort -V | tail -n 1)
    		CUDAFLAGS=$(pkg-config --cflags ${cuda_pkg})
    		CUDALDFLAGS=$(pkg-config --libs ${cuda_pkg})
            fi # end of CUDA runtime decision tree
            echo $CUDALDFLAGS | grep cublas > /dev/null
            cudart_check=$?
            if [ "x${cudart_check}" != "x0" ]; then
		if [ "x${CUBLAS}" != "xno" ]; then
		    CUDALDFLAGS="$CUDALDFLAGS -lcublas"
		fi
            fi
            echo $CUDALDFLAGS | grep cudart > /dev/null
            cudart_check=$?
            if [ "x${cudart_check}" != "x0" ]; then
		if [ "x${CUBLAS}" != "xno" ]; then
		    CUDALDFLAGS="$CUDALDFLAGS -lcudart"
		fi
            fi
	else
            # pkg-config is not available
            if [ -f /usr/local/cuda/include/cuda.h ]; then
		CUDAFLAGS="-I/usr/local/cuda/include"
		CUDALDFLAGS="-L/usr/local/cuda/lib64 -lcublas -lcudart"
            elif [ -f /usr/include/cuda.h ]; then
		CUDAFLAGS="-I/usr/include"
		CUDALDFLAGS="-lcublas -lcudart"
            elif [ "x$CUDA_HOME" != "x" ]; then
		CUDAFLAGS="-I${CUDA_HOME}/include"
		CUDALDFLAGS="-L${CUDA_HOME}/lib64 -lcublas -lcudart"
            fi
	fi # end of pkg-config decision tree
    fi # end of CUDAFLAGS user override protection
    if [ "x$CUDAFLAGS$CUDALDFLAGS" != "x" ]; then
	# WARNING: the corresponding test in configure.ac has an error!
	CUBLASFLAGS=" -DUSE_CUBLAS ${CUDAFLAGS}"
	CUBLASLDFLAGS="${CUDALDFLAGS}"
    fi
    if [ "x$CUBLASFLAGS$CUBLASLDFLAGS" != "x" ]; then
	# The matrices are inverted with cuSOLVER (part of the CUDA toolkit;
	# its 64-bit interface needs CUDA 11.1 or later): cuBLAS is used only
	# if cuSOLVER is available too.
	cat > test_cusolver.cpp <<EOF
#include <cusolverDn.h>
int main() {
  cusolverDnHandle_t handle = 0;
  cusolverDnParams_t params = 0;
  size_t device_bytes = 0, host_bytes = 0;
  if (cusolverDnCreate(&handle) != CUSOLVER_STATUS_SUCCESS) return 0;
  cusolverDnCreateParams(&params);
  cusolverDnXgetrf_bufferSize(handle, params, 2, 2, CUDA_C_64F, 0, 2, CUDA_C_64F, &device_bytes, &host_bytes);
  cusolverDnXgetrs(handle, params, CUBLAS_OP_N, 2, 2, CUDA_C_64F, 0, 2, 0, CUDA_C_64F, 0, 2, 0);
  cusolverDnDestroyParams(params);
  cusolverDnDestroy(handle);
  return 0;
}
EOF
	$CXX${CLANGFLAGS} ${CUDAFLAGS} test_cusolver.cpp -o test_cusolver ${CUDALDFLAGS} -lcusolver > /dev/null 2>>error.log
	if [ "x$?" = "x0" ]; then
	    CUBLASLDFLAGS="${CUBLASLDFLAGS} -lcusolver"
	else
	    CUBLASFLAGS=""
	    CUBLASLDFLAGS=""
	    CUBLAS_NO_CUSOLVER="yes"
	fi
	rm -f test_cusolver test_cusolver.cpp
    fi
    if [ "x$CUBLASFLAGS$CUBLASLDFLAGS" = "x" ]; then
	if [ "x$CUBLAS_NO_CUSOLVER" = "xyes" ]; then
	    echo "no (cuSOLVER not found)"
	    echo "no (cuSOLVER not found)" >>configure.log
	else
	    echo "no"
	    echo "no" >>configure.log
	fi
	if [ "x$CUBLAS" = "xyes" ]; then
	    echo "ERROR: cuBLAS was required, but no cuBLAS with cuSOLVER found."
	    echo "ERROR: cuBLAS was required, but no cuBLAS with cuSOLVER found." >>configure.log
	    exit 2
	fi
    else
	echo "cuBLAS, cuSOLVER"
	echo "cuBLAS, cuSOLVER" >>configure.log
    fi
else
    CUBLASFLAGS=""
    CUBLASLDFLAGS=""
fi
# End of cuBLAS checks
# MAGMA checks
if [ "x$MAGMA" != "xno" ]; then
    echo -n "configure: checking for MAGMA... "
    echo -n "configure: checking for MAGMA... " >>configure.log
    if [ "x$ENABLE_ILP64" = "xyes" ]; then
	# 64-bit indices are enabled
	MAGMA_ILP64_FLAG="-DMAGMA_ILP64"
	MAGMA_LD_SPEC="64"
    else
	# 64-bit indices are disabled
	MAGMA_ILP64_FLAG=""
	MAGMA_LD_SPEC=""
    fi # end of 64-bit decision tree
    pkg-config --version > /dev/null
    use_pkg_config=$?
    # Find the MAGMA headers and library first: the GPU runtime flags depend
    # on the back end MAGMA was built on, CUDA or HIP.
    MAGMA_INCLUDE=""
    MAGMA_SOURCE="none"
    if [ "x${MAGMA_ROOT}${MAGMA_HOME}${MAGMA_DIR}" = "x" ]; then
	# MAGMA environment is not defined
	if [ "x$use_pkg_config" = "x0" ]; then
            # use pkg-config to search for MAGMA
            declare -a pkg_array=$(pkg-config --list-all | grep magma)
            for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep magma > /dev/null
            result=$?
	    if [ "x$result" = "x0" ]; then
		# MAGMA was found
		magma_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep magma)
		MAGMA_INCLUDE=$(pkg-config --cflags-only-I ${magma_pkg})
		MAGMA_LIBS_DIR=$(pkg-config --libs-only-L ${magma_pkg})
		MAGMALDFLAGS=" ${MAGMA_LIBS_DIR} -lmagma"
		MAGMA_SOURCE="pkg-config"
	    fi # end of MAGMA decision tree
	else
            # search for MAGMA in some standard folders
            if [ -f /usr/include/magma_v2.h ]; then
		MAGMA_INCLUDE="-I/usr/include"
		MAGMALDFLAGS=" -lmagma"
		MAGMA_SOURCE="folders"
            elif [ -f /usr/local/include/magma_v2.h ]; then
		MAGMA_INCLUDE="-I/usr/local/include"
		MAGMALDFLAGS=" -lmagma"
		MAGMA_SOURCE="folders"
            fi
	fi # end of pkg-config decision tree
    else
	# MAGMA environment is defined, so configuration makes sure that MAGMA_ROOT is defined too
	if [ "x${MAGMA_HOME}" != "x" ]; then
            MAGMA_ROOT="${MAGMA_HOME}"
	elif [ "x${MAGMA_DIR}" != "x" ]; then
            MAGMA_ROOT="${MAGMA_DIR}"
	fi
	MAGMA_INCLUDE="-I${MAGMA_ROOT}/include"
	MAGMALDFLAGS=" -L${MAGMA_ROOT}/lib -lmagma"
	MAGMA_SOURCE="root"
    fi
    # The back end: magma_config.h, generated when MAGMA is built, defines
    # MAGMA_HAVE_CUDA or MAGMA_HAVE_HIP. Without it, assume CUDA.
    GPURTFLAGS=""
    GPURTLDFLAGS=""
    if [ "x$MAGMA_SOURCE" != "xnone" ]; then
	MAGMA_BACKEND="cuda"
	for magma_dir in $(echo "$MAGMA_INCLUDE" | tr ' ' '\n' | sed -n 's/^-I//p'); do
	    if [ -f "$magma_dir/magma_config.h" ]; then
		grep -q "^[[:space:]]*#[[:space:]]*define[[:space:]][[:space:]]*MAGMA_HAVE_HIP" "$magma_dir/magma_config.h"
		if [ "x$?" = "x0" ]; then
		    MAGMA_BACKEND="hip"
		fi
		break
	    fi
	done
    fi
    if [ "x$MAGMA_BACKEND" = "xhip" ]; then
	# MAGMA on HIP (AMD GPUs): the HIP runtime, from ROCm.
	rocm_root="${ROCM_PATH:-/opt/rocm}"
	if [ -f "${rocm_root}/include/hip/hip_runtime_api.h" ]; then
	    GPURTFLAGS="-D__HIP_PLATFORM_AMD__ -I${rocm_root}/include"
	    GPURTLDFLAGS="-L${rocm_root}/lib -lamdhip64"
	else
	    echo -n "(MAGMA uses HIP, but no HIP headers in ${rocm_root}: set ROCM_PATH) "
	    echo -n "(MAGMA uses HIP, but no HIP headers in ${rocm_root}: set ROCM_PATH) " >>configure.log
	    MAGMA_SOURCE="none"
	fi
	if [ "x$CUBLASFLAGS$CUBLASLDFLAGS" != "x" ]; then
	    # cuBLAS is CUDA: it cannot share the program with a MAGMA on HIP.
	    if [ "x$CUBLAS" = "xyes" ]; then
		echo ""
		echo "ERROR: cuBLAS was required, but MAGMA was built on HIP (AMD GPUs)."
		echo "ERROR: cuBLAS was required, but MAGMA was built on HIP (AMD GPUs)." >>configure.log
		exit 2
	    fi
	    CUBLASFLAGS=""
	    CUBLASLDFLAGS=""
	    CUBLAS_DROPPED="yes"
	fi
    elif [ "x$MAGMA_BACKEND" = "xcuda" ]; then
	# MAGMA on CUDA (NVIDIA GPUs): the CUDA headers and runtime library.
	if [ "x${CUDAFLAGS}" = "x" ]; then
	    # CUDA_HOME or CUDA_PATH select the CUDA installation; otherwise the
	    # newest pkg-config package, otherwise the standard folders.
	    if [ "x$CUDA_HOME" != "x" ]; then
		CUDAFLAGS="-I${CUDA_HOME}/include"
	    elif [ "x$CUDA_PATH" != "x" ]; then
		CUDAFLAGS="-I${CUDA_PATH}/include"
	    elif [ "x$use_pkg_config" = "x0" ]; then
		declare -a pkg_array=$(pkg-config --list-all | grep cudart)
		for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart > /dev/null
		result=$?
		if [ "x$result" = "x0" ]; then
		    # CUDA runtime detected
		    cuda_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart | sort -V | tail -n 1)
		    CUDAFLAGS=$(pkg-config --cflags ${cuda_pkg})
		    GPURTLDFLAGS=$(pkg-config --libs ${cuda_pkg})
		fi # end of CUDA runtime decision tree
	    fi
	    if [ "x$CUDAFLAGS" = "x" ]; then
		if [ -f /usr/local/cuda/include/cuda.h ]; then
		    CUDAFLAGS="-I/usr/local/cuda/include"
		elif [ -f /usr/include/cuda.h ]; then
		    CUDAFLAGS="-I/usr/include"
		fi
	    fi
	fi # end of CUDAFLAGS user override protection
	if [ "x$CUDAFLAGS" = "x" -a "x$MAGMA_SOURCE" = "xfolders" ]; then
	    # A MAGMA in the standard folders is used only with the CUDA headers.
	    MAGMA_SOURCE="none"
	fi
	if [ "x$GPURTLDFLAGS" = "x" ]; then
	    # The programs call the CUDA runtime (gpu_runtime.cpp) directly.
	    cuda_rt_root=$(echo "$CUDAFLAGS" | sed -n 's/.*-I\([^ ]*\)\/include.*/\1/p')
	    if [ "x$cuda_rt_root" != "x" -a -d "${cuda_rt_root}/lib64" ]; then
		GPURTLDFLAGS="-L${cuda_rt_root}/lib64 -lcudart"
	    elif [ "x$cuda_rt_root" != "x" -a -d "${cuda_rt_root}/lib" ]; then
		GPURTLDFLAGS="-L${cuda_rt_root}/lib -lcudart"
	    else
		GPURTLDFLAGS="-lcudart"
	    fi
	fi
	GPURTFLAGS="$CUDAFLAGS"
    fi
    if [ "x$MAGMA_SOURCE" != "xnone" ]; then
	MAGMAFLAGS=" -DUSE_MAGMA ${MAGMA_ILP64_FLAG} ${GPURTFLAGS} ${MAGMA_INCLUDE}"
	MAGMALDFLAGS="${MAGMALDFLAGS} ${GPURTLDFLAGS}"
    else
	MAGMA_BACKEND="none"
	MAGMAFLAGS=""
	MAGMALDFLAGS=""
    fi
    if [ "x$MAGMAFLAGS$MAGMALDFLAGS" = "x" ]; then
	echo "no"
	echo "no" >>configure.log
	if [ "x$MAGMA" = "xyes" ]; then
	    echo "ERROR: MAGMA was required, but no MAGMA found."
	    echo "ERROR: MAGMA was required, but no MAGMA found." >>configure.log
	    exit 2
	fi
    else
	if [ "x$MAGMAFLAGS" = "x" ]; then
	    MAGMAFLAGS="$MAGMA_INVERT_FLAGS"
	else
	    MAGMAFLAGS="$MAGMAFLAGS$MAGMA_INVERT_FLAGS"
	fi
	if [ "x$MAGMA_BACKEND" = "xhip" ]; then
	    echo "yes (HIP)"
	    echo "yes (HIP)" >>configure.log
	else
	    echo "yes (CUDA)"
	    echo "yes (CUDA)" >>configure.log
	fi
	if [ "x$CUBLAS_DROPPED" = "xyes" ]; then
	    echo "WARNING: cuBLAS is not used, because MAGMA was built on HIP (AMD GPUs)."
	    echo "WARNING: cuBLAS is not used, because MAGMA was built on HIP (AMD GPUs)." >>configure.log
	fi
	# The size of magma_int_t is fixed when MAGMA is built, while magma_v2.h
	# follows MAGMA_ILP64: check it on the library, by letting MAGMA write
	# integers whose upper 32 bits are preset.
	echo -n "configure: checking whether MAGMA uses 64-bit integers... "
	echo -n "configure: checking whether MAGMA uses 64-bit integers... " >>configure.log
	cat > test_magma.cpp <<EOF
#include <magma_v2.h>
int main() {
  magma_int_t major = ((magma_int_t)7) << 32;
  magma_int_t minor = major, micro = major;
  magma_version(&major, &minor, &micro);
  // 0: MAGMA wrote 64-bit integers; 2: it wrote only the lower 32 bits.
  return ((major >> 32) == 7) ? 2 : 0;
}
EOF
	magma_ints="unknown"
	$CXX${CLANGFLAGS} ${MAGMAFLAGS} -DMAGMA_ILP64 test_magma.cpp -o test_magma ${MAGMALDFLAGS} ${CUBLASLDFLAGS} > /dev/null 2>>error.log
	if [ "x$?" = "x0" ]; then
	    ./test_magma > /dev/null 2>>error.log
	    result=$?
	    if [ "x$result" = "x0" ]; then
		magma_ints="64"
	    elif [ "x$result" = "x2" ]; then
		magma_ints="32"
	    fi
	fi
	rm -f test_magma test_magma.cpp
	if [ "x$magma_ints" = "x64" ]; then
	    echo "yes"
	    echo "yes" >>configure.log
	    MAGMAFLAGS=$(echo "$MAGMAFLAGS" | sed 's/ -DMAGMA_ILP64//g')
	    MAGMAFLAGS=$(echo "$MAGMAFLAGS" | sed 's/-DUSE_MAGMA/-DUSE_MAGMA -DMAGMA_ILP64/')
	elif [ "x$magma_ints" = "x32" ]; then
	    echo "no"
	    echo "no" >>configure.log
	    MAGMAFLAGS=$(echo "$MAGMAFLAGS" | sed 's/ -DMAGMA_ILP64//g')
	    echo "$LAPACKFLAGS" | grep -q "LAPACK_ILP64"
	    if [ "x$?" = "x0" ]; then
		echo "WARNING: MAGMA uses 32-bit integers, while LAPACK is used with 64-bit integers."
		echo "WARNING: MAGMA uses 32-bit integers, while LAPACK is used with 64-bit integers." >>configure.log
	    fi
	else
	    echo "unknown (keeping the 64-bit integer setting)"
	    echo "unknown (keeping the 64-bit integer setting)" >>configure.log
	fi
    fi
else
    MAGMAFLAGS=""
    MAGMALDFLAGS=""
fi
# End of MAGMA checks
# Offload checks
if [ "x$OFFLOAD" != "xno" ]; then
    # Offloading is configured differently for clang and for g++, and for
    # NVIDIA (nvptx) and AMD (amdgcn) GPUs. clang needs the target triple, the
    # GPU architecture and a CUDA installation (ptxas, libdevice) or a ROCm
    # one; g++ needs its nvptx-none or amdgcn-amdhsa offload compiler.
    OFFLOAD_COMPILER="gcc"
    $CXX --version | grep "clang" > /dev/null
    if [ "x$?" = "x0" ]; then
	OFFLOAD_COMPILER="clang"
    fi
    amd_gpu_arch=""
    if [ "x$OFFLOAD_TARGET" != "xnvptx" ]; then
	amd_gpu_arch=$(amdgpu-arch 2>>error.log | head -n 1)
	if [ "x$amd_gpu_arch" = "x" ]; then
	    amd_gpu_arch=$(rocminfo 2>>error.log | grep -o "gfx[0-9a-f]*" | head -n 1)
	fi
    fi
    nv_gpu_cap=""
    if [ "x$OFFLOAD_TARGET" != "xamdgcn" ]; then
	nv_gpu_cap=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>>error.log | head -n 1 | tr -d ' .')
    fi
    if [ "x$OFFLOAD_TARGET" = "xauto" ]; then
	# The GPUs of the GPU libraries come first, then the GPUs found here.
	if [ "x$MAGMA_BACKEND" = "xhip" ]; then
	    OFFLOAD_TARGET="amdgcn"
	elif [ "x$MAGMA_BACKEND" = "xcuda" -o "x$CUBLASFLAGS$CUBLASLDFLAGS" != "x" ]; then
	    OFFLOAD_TARGET="nvptx"
	elif [ "x$nv_gpu_cap" != "x" ]; then
	    OFFLOAD_TARGET="nvptx"
	elif [ "x$amd_gpu_arch" != "x" ]; then
	    OFFLOAD_TARGET="amdgcn"
	else
	    OFFLOAD_TARGET="nvptx"
	fi
    fi
    if [ "x$OFFLOAD_TARGET" = "xamdgcn" -a "x$MAGMA_BACKEND" = "xcuda" ]; then
	echo "WARNING: offloading to AMD GPUs, while MAGMA was built on CUDA (NVIDIA GPUs)."
	echo "WARNING: offloading to AMD GPUs, while MAGMA was built on CUDA (NVIDIA GPUs)." >>configure.log
    elif [ "x$OFFLOAD_TARGET" = "xnvptx" -a "x$MAGMA_BACKEND" = "xhip" ]; then
	echo "WARNING: offloading to NVIDIA GPUs, while MAGMA was built on HIP (AMD GPUs)."
	echo "WARNING: offloading to NVIDIA GPUs, while MAGMA was built on HIP (AMD GPUs)." >>configure.log
    fi
    if [ "x$OFFLOAD_ARCH" = "xauto" ]; then
	echo -n "configure: checking GPU architecture... "
	echo -n "configure: checking GPU architecture... " >>configure.log
	if [ "x$OFFLOAD_TARGET" = "xamdgcn" -a "x$amd_gpu_arch" != "x" ]; then
	    OFFLOAD_ARCH="$amd_gpu_arch"
	elif [ "x$OFFLOAD_TARGET" = "xnvptx" -a "x$nv_gpu_cap" != "x" -a "x$OFFLOAD_COMPILER" = "xclang" ]; then
	    OFFLOAD_ARCH="sm_${nv_gpu_cap}"
	fi
	if [ "x$OFFLOAD_ARCH" != "xauto" ]; then
	    echo "$OFFLOAD_ARCH"
	    echo "$OFFLOAD_ARCH" >>configure.log
	else
	    echo "unknown"
	    echo "unknown" >>configure.log
	    if [ "x$OFFLOAD_COMPILER" = "xclang" ]; then
		# Let clang query the GPU itself.
		OFFLOAD_ARCH="native"
	    elif [ "x$OFFLOAD_TARGET" = "xamdgcn" ]; then
		OFFLOAD_ARCH="gfx90a"
	    else
		# The architectures known to the nvptx back end depend on the
		# g++ version: keep the long-standing default.
		OFFLOAD_ARCH="sm_70"
	    fi
	fi
    fi
    cat > conf_test_offload.cpp <<EOF
#include <omp.h>

int main(int argc, char** argv) {
  const int size = 1000000;
  int *numbers = new int[size]();
  int on_host = 1;
#pragma omp target map(from:on_host)
  on_host = omp_is_initial_device();
#pragma omp target teams distribute parallel for map(tofrom:numbers[0:size])
  for (int i = 0; i < size; i++) {
    numbers[i] = 1;
  }
  int result = (numbers[0] == 1 && numbers[size - 1] == 1) ? 0 : 1;
  delete[] numbers;
  // The build machine may have no GPU (e.g. a cluster login node): the
  // test passes if the target region runs, on the GPU or on the host.
  return result + (on_host ? 2 : 0);
}
EOF
    result=1
    OFFLOAD_ON_HOST="no"
    if [ "x$OFFLOAD_COMPILER" = "xclang" ]; then
	if [ "x$OFFLOAD_TARGET" = "xamdgcn" ]; then
	    OFFLOAD_CUDA_PATH=""
	elif [ "x$OFFLOAD_CUDA_PATH" = "xauto" ]; then
	    # Look for a complete CUDA installation: clang needs ptxas and the
	    # libdevice bitcode library. The installation that gave the CUDA
	    # flags is preferred, so that offloading and the CUDA libraries use
	    # the same version.
	    OFFLOAD_CUDA_PATH=""
	    cuda_flags_root=$(echo "$CUDAFLAGS" | sed -n 's/.*-I\([^ ]*\)\/targets\/.*/\1/p')
	    if [ "x$cuda_flags_root" = "x" ]; then
		cuda_flags_root=$(echo "$CUDAFLAGS" | sed -n 's/.*-I\([^ ]*\)\/include.*/\1/p')
	    fi
	    cuda_candidates="${CUDA_HOME} ${CUDA_PATH} ${cuda_flags_root}"
	    nvcc_path=$(which nvcc 2>/dev/null)
	    if [ "x$nvcc_path" != "x" ]; then
		cuda_candidates="$cuda_candidates $(dirname $(dirname $(readlink -f $nvcc_path)))"
	    fi
	    cuda_candidates="$cuda_candidates /usr/local/cuda $(ls -d /usr/local/cuda-* 2>/dev/null | sort -V -r)"
	    for cuda_dir in $cuda_candidates; do
		if [ -x "$cuda_dir/bin/ptxas" -a -d "$cuda_dir/nvvm/libdevice" ]; then
		    OFFLOAD_CUDA_PATH="$cuda_dir"
		    break
		fi
	    done
	fi
	CUDAPATHFLAGS=""
	if [ "x$OFFLOAD_CUDA_PATH" != "x" ]; then
	    CUDAPATHFLAGS=" --cuda-path=${OFFLOAD_CUDA_PATH}"
	fi
	if [ "x$OFFLOAD_TARGET" = "xamdgcn" ]; then
	    offload_triple="amdgcn-amd-amdhsa"
	    if [ -d "${ROCM_PATH:-/opt/rocm}" ]; then
		CUDAPATHFLAGS=" --rocm-path=${ROCM_PATH:-/opt/rocm}"
	    fi
	else
	    offload_triple="nvptx64-nvidia-cuda"
	fi
	echo -n "configure: checking whether $CXX supports offload to ${OFFLOAD_ARCH}... "
	echo -n "configure: checking whether $CXX supports offload to ${OFFLOAD_ARCH}... " >>configure.log
	# Two spellings of the target options, for different clang versions.
	for target_flags in \
	    " -fopenmp-targets=${offload_triple} --offload-arch=${OFFLOAD_ARCH}" \
	    " -fopenmp-targets=${offload_triple} -Xopenmp-target=${offload_triple} -march=${OFFLOAD_ARCH}"; do
	    echo "configure: trying $CXX -fopenmp${target_flags}${CUDAPATHFLAGS}" >>error.log
	    $CXX${CLANGFLAGS} -fopenmp${target_flags}${CUDAPATHFLAGS} conf_test_offload.cpp -o conf_test_offload > /dev/null 2>>error.log
	    result=$?
	    if [ "x$result" = "x0" ]; then
		./conf_test_offload > /dev/null 2>>error.log
		result=$?
		echo "configure: offload test program exited with $result" >>error.log
		rm -f conf_test_offload
		if [ "x$result" = "x2" ]; then
		    OFFLOAD_ON_HOST="yes"
		    result=0
		fi
	    fi
	    if [ "x$result" = "x0" ]; then
		OFFLOADFLAGS=" -DUSE_TARGET_OFFLOAD${target_flags}${CUDAPATHFLAGS}"
		# clang uses its own OpenMP runtime: do not link libgomp.
		OFFLOADLDFLAGS=""
		break
	    fi
	done
    else
	echo -n "configure: checking whether $CXX supports offload to ${OFFLOAD_ARCH}... "
	echo -n "configure: checking whether $CXX supports offload to ${OFFLOAD_ARCH}... " >>configure.log
	if [ "x$OFFLOAD_TARGET" = "xamdgcn" ]; then
	    gcc_offload_target="amdgcn-amdhsa"
	    gcc_offload_options="-O2 -march=${OFFLOAD_ARCH}"
	else
	    gcc_offload_target="nvptx-none"
	    gcc_offload_options="-O2 -march=${OFFLOAD_ARCH} -mptx=7.3"
	fi
	$CXX -fopenmp -fcf-protection=none -fno-stack-protector -fno-strict-aliasing -foffload=${gcc_offload_target}="${gcc_offload_options}" conf_test_offload.cpp -o conf_test_offload > /dev/null 2>>error.log
	result=$?
	if [ "x$result" = "x0" ]; then
	    ./conf_test_offload > /dev/null 2>>error.log
	    result=$?
	    echo "configure: offload test program exited with $result" >>error.log
	    rm -f conf_test_offload
	    if [ "x$result" = "x2" ]; then
		OFFLOAD_ON_HOST="yes"
		result=0
	    fi
	fi
	if [ "x$result" = "x0" ]; then
	    OFFLOADFLAGS=" -DUSE_TARGET_OFFLOAD -fcf-protection=none -fno-stack-protector -fno-strict-aliasing -foffload=${gcc_offload_target}=\"${gcc_offload_options} -lm\""
	    OFFLOADLDFLAGS=" -lgomp"
	fi
    fi
    rm -f conf_test_offload.cpp
    if [ "x$result" = "x0" ]; then
	echo "yes"
	echo "yes" >>configure.log
	if [ "x$OFFLOAD_ON_HOST" = "xyes" ]; then
	    echo "WARNING: the offload test ran on the host: no GPU was available to it."
	    echo "WARNING: the offload test ran on the host: no GPU was available to it." >>configure.log
	fi
	if [ "x${OMPFLAGS}" = "x" ]; then
	    OFFLOADFLAGS=" -fopenmp${OFFLOADFLAGS}"
	fi
    else
	echo "no"
	echo "no" >>configure.log
	OFFLOADFLAGS=""
	OFFLOADLDFLAGS=""
    fi
else
    OFFLOADFLAGS=""
    OFFLOADLDFLAGS=""
fi
# End of offload checks
if [ "x$CXXFLAGS" = "x" ]; then
    CXXFLAGS="-O${CXX_OPT}${CXX_DBG}${CLANGFLAGS}${INCLUDEFLAGS}${HDF5FLAGS}${OMPFLAGS}${MPIFLAGS}${LAPACKFLAGS}${CUBLASFLAGS}${MAGMAFLAGS}${DEBUGFLAGS}${OFFLOADFLAGS}${NVTXFLAGS}"
fi
# The version is added also to user-defined CXXFLAGS.
CXXFLAGS="${CXXFLAGS}${VERSIONFLAGS}"
if [ "x$CXXLDFLAGS" = "x" ]; then
    if [ "x$LIBMODE" = "xstatic" ]; then
	CXXLDFLAGS="${ASNEEDED}-Llibnptm -lnptm ${HDF5LDFLAGS} ${LDFLAGS} ${LAPACKLDFLAGS} ${CUBLASLDFLAGS} ${MAGMALDFLAGS} ${OFFLOADLDFLAGS}${OMPLDFLAGS}"
    else
	CXXLDFLAGS="${ASNEEDED}-Llibnptm -lnptm ${HDF5LDFLAGS} ${LDFLAGS} ${LAPACKLDFLAGS} ${CUBLASLDFLAGS} ${MAGMALDFLAGS} ${OFFLOADLDFLAGS}${OMPLDFLAGS}"
    fi
fi

if [ "x$FFLAGS" != "x" ]; then
    FCFLAGS=$FFLAGS
fi
# End of configuration logic

# Print a summary of configuration options
echo "INFO: optimization level is ${CXX_OPT}."
echo "INFO: optimization level is ${CXX_OPT}." >>configure.log
if [ "x${CXX_DBG}" = "x" ]; then
    echo "INFO: gdb is disabled."
    echo "INFO: gdb is disabled." >>configure.log
else
    echo "INFO: gdb is enabled."
    echo "INFO: gdb is enabled." >>configure.log
fi
if [ "x${OMPFLAGS}" = "x" ]; then
    echo "INFO: OpenMP is disabled."
    echo "INFO: OpenMP is disabled." >>configure.log
else
    echo "INFO: OpenMP is enabled."
    echo "INFO: OpenMP is enabled." >>configure.log
fi
if [ "x${MPIFLAGS}" = "x" ]; then
    echo "INFO: MPI is disabled."
    echo "INFO: MPI is disabled." >>configure.log
else
    echo "INFO: MPI is enabled."
    echo "INFO: MPI is enabled." >>configure.log
fi
if [ "x${LAPACKFLAGS}" = "x" ]; then
    echo "INFO: LAPACK is disabled."
    echo "INFO: LAPACK is disabled." >>configure.log
else
    echo "INFO: LAPACK is enabled."
    echo "INFO: LAPACK is enabled." >>configure.log
fi
if [ "x${CUBLASFLAGS}" = "x" ]; then
    echo "INFO: cuBLAS was not found."
    echo "INFO: cuBLAS was not found." >>configure.log
else
    echo "INFO: cuBLAS was found (inversion with cuSOLVER)."
    echo "INFO: cuBLAS was found (inversion with cuSOLVER)." >>configure.log
fi
if [ "x${MAGMAFLAGS}" = "x" ]; then
    echo "INFO: MAGMA is disabled."
    echo "INFO: MAGMA is disabled." >>configure.log
else
    echo "INFO: MAGMA is enabled."
    echo "INFO: MAGMA is enabled." >>configure.log
    if [ "x${MAGMA_INVERT_FLAGS}" = "x" ]; then
	echo "INFO: using LU factorisation for matrix inversion."
	echo "INFO: using LU factorisation for matrix inversion." >>configure.log
    elif [ "x${MAGMA_INVERT_FLAGS}" = "x -DUSE_ZGESV_GPU" ]; then
	echo "INFO: using MAGMA zgesv_gpu function for matrix inversion."
	echo "INFO: using MAGMA zgesv_gpu function for matrix inversion." >>configure.log
    elif [ "x${MAGMA_INVERT_FLAGS}" = "x -DUSE_ZGESV_RBT" ]; then
	echo "INFO: using MAGMA zgesv_rbt function for matrix inversion."
	echo "INFO: using MAGMA zgesv_rbt function for matrix inversion." >>configure.log
    fi
fi
if [ "x${NVTXFLAGS}" = "x" ]; then
    echo "INFO: NVTX profiling is disabled."
    echo "INFO: NVTX profiling is disabled." >>configure.log
else
    echo "INFO: NVTX profiling is enabled."
    echo "INFO: NVTX profiling is enabled." >>configure.log
fi
if [ "x$MAGMA_BACKEND" = "xhip" -a "x$MAGMAFLAGS" != "x" ] || [ "x${OFFLOADFLAGS}" != "x" -a "x$OFFLOAD_TARGET" = "xamdgcn" ]; then
    echo "WARNING: support for AMD GPUs (MAGMA on HIP, offload to amdgcn) is experimental and not yet tested."
    echo "WARNING: support for AMD GPUs (MAGMA on HIP, offload to amdgcn) is experimental and not yet tested." >>configure.log
fi
if [ "x${OFFLOADFLAGS}" = "x" ]; then
    echo "INFO: GPU offload through OpenMP is disabled."
    echo "INFO: GPU offload through OpenMP is disabled." >>configure.log
else
    echo "INFO: GPU offload through OpenMP is enabled (${OFFLOAD_COMPILER}, ${OFFLOAD_TARGET}, ${OFFLOAD_ARCH})."
    echo "INFO: GPU offload through OpenMP is enabled (${OFFLOAD_COMPILER}, ${OFFLOAD_TARGET}, ${OFFLOAD_ARCH})." >>configure.log
fi
if [ "x${LIBMODE}" = "xstatic" ]; then
    echo "INFO: configured to build static proprietary libraries."
    echo "INFO: configured to build static proprietary libraries." >>configure.log
else
    echo "INFO: configured to build shared proprietary libraries."
    echo "INFO: configured to build shared proprietary libraries." >>configure.log
fi
# End of summary printing section.

# Write all the configuration settings to Makefile.in
cat > Makefile.in <<EOF
AR=${AR}
FC=${FC}
FCFLAGS=${FCFLAGS}
BUILDFORTRAN=${BUILDFORTRAN}
LIBMODE=${LIBMODE}
CXX=${CXX}
CXXFLAGS=${CXXFLAGS}
CXXLDFLAGS=${CXXLDFLAGS}
EOF
# The compiler behind an MPI wrapper is chosen by the environment; record it,
# so that make uses the compiler that configure tested.
if [ "x${OMPI_CXX}" != "x" ]; then
    echo "export OMPI_CXX=${OMPI_CXX}" >> Makefile.in
fi
if [ "x${MPICH_CXX}" != "x" ]; then
    echo "export MPICH_CXX=${MPICH_CXX}" >> Makefile.in
fi

echo "configure: finished."
echo "configure: finished." >>configure.log
# End of script execution
