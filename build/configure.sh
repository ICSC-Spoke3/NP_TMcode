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
    echo "--enable-openmp                   Enable OpenMP multi-threading (DEFAULT).       "
    echo "--disable-openmp                  Disable OpenMP multi-threading.                "
    echo "--enable-debug=FEATURE            Enable debug output of specified feature.      "
    echo "--enable-nvtx                     Enable NVTX profiling tools.                   "
    echo "--disable-nvtx                    Disable NVTX profiling tools (DEFAULT).        "
    echo "--enable-optimize=OPT             Use OPT-level compiler optimization.           "
    echo "--enable-shared                   Use shared libraries (default is static).      "
    echo "--disable-shared                  Use static libraries (DEFAULT).                "
    echo "--help                            Print this help and exit.                      "
    echo "--with-cublas                     Use cuBLAS (DEFAULT).                          "
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
    echo "CUDAFLAGS                         CUDA compilation flags.                        "
    echo "CUDALDFLAGS                       CUDA libraries linking flags.                  "
    echo "CXX                               C++ compiler.                                  "
    echo "CXXFLAGS                          C++ compilation flags.                         "
    echo "CXXLDFLAGS                        C++ linker flags.                              "
    echo "FC                                FORTRAN compiler.                              "
    echo "HDF5_INCLUDE                      Path to the HDF5 header files.                 "
    echo "HDF5_LIB                          Path to the HDF5 library files.                "
    echo "MAGMA_HOME                        Root folder of MAGMA distribution.             "
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
            LAPACKFLAGS=" -DUSE_LAPACK -DUSE_MKL ${LAPACK_ILP64_FLAG} -I{MKLROOT}/include"
            LAPACKLDFLAGS=" -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel${LAPACK_ILP64_LDSPEC} -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
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
    if [ "x${CUDAFLAGS}${CUDALDFLAGS}" = "x" ]; then
	if [ "x$use_pkg_config" = "x0" ]; then
            # pkg-config is available
            declare -a pkg_array=$(pkg-config --list-all | grep cublas)
            for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cublas > /dev/null
            result=$?
            if [ "x$result" = "x0" ]; then
		# CUBLAS detected
		cuda_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cublas)
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
		cuda_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart)
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
    if [ "x$CUBLASFLAGS$CUBLASLDFLAGS" = "x" ]; then
	echo "no"
	echo "no" >>configure.log
	if [ "x$CUBLAS" = "xyes" ]; then
	    echo "ERROR: cuBLAS was required, but no cuBLAS found."
	    echo "ERROR: cuBLAS was required, but no cuBLAS found." >>configure.log
	    exit 2
	fi
    else
	echo "cuBLAS"
	echo "cuBLAS" >>configure.log
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
    if [ "x${CUDAFLAGS}" = "x" ]; then
	if [ "x$use_pkg_config" = "x0" ]; then
            # pkg-config is available
            declare -a pkg_array=$(pkg-config --list-all | grep cudart)
            for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart > /dev/null
            result=$?
            if [ "x$result" = "x0" ]; then
		# CUDA runtime detected
		cuda_pkg=$(for i in "${pkg_array[@]}"; do echo "$i" | cut --delimiter=" " -f1; done | grep cudart)
    		CUDAFLAGS=$(pkg-config --cflags ${cuda_pkg})
            fi # end of CUDA runtime decision tree
	else
            # pkg-config is not available
            if [ -f /usr/local/cuda/include/cuda.h ]; then
		CUDAFLAGS="-I/usr/local/cuda/include"
            elif [ -f /usr/include/cuda.h ]; then
		CUDAFLAGS="-I/usr/include"
            elif [ "x$CUDA_HOME" != "x" ]; then
		CUDAFLAGS="-I${CUDA_HOME}/include"
            fi
	fi # end of pkg-config decision tree
    fi # end of CUDAFLAGS user override protection
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
		MAGMAFLAGS=" -DUSE_MAGMA ${MAGMA_ILP64_FLAG} $CUDAFLAGS ${MAGMA_INCLUDE}"
		MAGMALDFLAGS=" ${MAGMA_LIBS_DIR}-lmagma"
	    fi # end of MAGMA decision tree
	else
            # search for MAGMA in some standard folders
        if [ "x$CUDAFLAGS" != "x" ]; then
            if [ -f /usr/include/magma_v2.h ]; then
		MAGMAFLAGS=" -DUSE_MAGMA ${MAGMA_ILP64_FLAG} $CUDAFLAGS -I/usr/include"
		MAGMALDFLAGS=" -lmagma"
            elif [ -f /usr/local/include/magma_v2.h ]; then
		MAGMAFLAGS=" -DUSE_MAGMA ${MAGMA_ILP64_FLAG} $CUDAFLAGS -I/usr/local/include"
		MAGMALDFLAGS=" -lmagma"
            fi
        fi
	fi # end of pkg-config decision tree
    else
	# MAGMA environment is defined, so configuration makes sure that MAGMA_ROOT is defined too
	if [ "x${MAGMA_HOME}" != "x" ]; then
            MAGMA_ROOT="${MAGMA_HOME}"
	elif [ "x${MAGMA_DIR}" != "x" ]; then
            MAGMA_ROOT="${MAGMA_DIR}"
	fi
	MAGMAFLAGS=" -DUSE_MAGMA -DMAGMA_ILP64 $CUDAFLAGS -I${MAGMA_ROOT}/include"
	MAGMALDFLAGS=" -L${MAGMA_ROOT}/lib -lmagma"
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
	echo "yes"
	echo "yes" >>configure.log
    fi
else
    MAGMAFLAGS=""
    MAGMALDFLAGS=""
fi
# End of MAGMA checks
# Offload checks
if [ "x$OFFLOAD" != "xno" ]; then
    echo -n "configure: checking whether system supports offload... "
    echo -n "configure: checking whether system supports offload... " >>configure.log
    cat > conf_test_offload.cpp <<EOF
#include <omp.h>

#pragma omp begin declare target device_type(any)
void fill_with_ones(int *array) {
#pragma omp teams distribute parallel for collapse(2)
  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < 1000; j++) {
      array[(1000 * i) + j] = 1;
    }
  }
}
#pragma omp end declare target

int main(int argc, char** argv) {
  int *numbers = new int[1000000]();
#pragma omp target map(tofrom:numbers[0:1000000])
  fill_with_ones(numbers);
  delete[] numbers;
  return 0;
}
EOF
    $CXX -fopenmp -fcf-protection=none -fno-stack-protector -foffload=nvptx-none="-O3 -ggdb -fcf-protection=none -fno-stack-protector -fopt-info -lm -latomic -lgomp" conf_test_offload.cpp -o conf_test_offload > /dev/null 2>>error.log
    result=$?
    rm conf_test_offload.cpp
    if [ "x$result" = "x0" ]; then
	./conf_test_offload > /dev/null 2>>error.log
	result=$?
	rm conf_test_offload
    fi
    if [ "x$result" = "x0" ]; then
	echo "yes"
	echo "yes" >>configure.log
	OFFLOADFLAGS=" -DUSE_TARGET_OFFLOAD -fcf-protection=none -fno-stack-protector -foffload=nvptx-none=\"-O${CXX_OPT}${CXX_DBG} -fcf-protection=none -fno-stack-protector -fopt-info -lm -latomic -lgomp\""
	if [ "x${OMPFLAGS}" = "x" ]; then
	    OFFLOADFLAGS="-fopenmp ${OFFLOADFLAGS}"
	fi
    else
	echo "no"
	echo "no" >>configure.log
	OFFLOADFLAGS=""
    fi
else
    OFFLOADFLAGS=""
fi
# End of offload checks
if [ "x$CXXFLAGS" = "x" ]; then
    CXXFLAGS="-O${CXX_OPT}${CXX_DBG}${CLANGFLAGS}${INCLUDEFLAGS}${HDF5FLAGS}${OMPFLAGS}${MPIFLAGS}${LAPACKFLAGS}${CUBLASFLAGS}${MAGMAFLAGS}${DEBUGFLAGS}${OFFLOADFLAGS}${NVTXFLAGS}"
fi
if [ "x$CXXLDFLAGS" = "x" ]; then
    if [ "x$LIBMODE" = "xstatic" ]; then
	CXXLDFLAGS="-Llibnptm -lnptm ${HDF5LDFLAGS} ${LDFLAGS} ${LAPACKLDFLAGS}${CUBLASLDFLAGS}${MAGMALDFLAGS}"
    else
	CXXLDFLAGS="-Llibnptm -lnptm ${HDF5LDFLAGS} ${LDFLAGS} ${LAPACKLDFLAGS}${CUBLASLDFLAGS}${MAGMALDFLAGS}"
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
    echo "INFO: cuBLAS was found."
    echo "INFO: cuBLAS was found." >>configure.log
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
if [ "x${OFFLOADFLAGS}" = "x" ]; then
    echo "INFO: GPU offload through OpenMP is disabled."
    echo "INFO: GPU offload through OpenMP is disabled." >>configure.log
else
    echo "INFO: GPU offload through OpenMP is enabled."
    echo "INFO: GPU offload through OpenMP is enabled." >>configure.log
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

echo "configure: finished."
echo "configure: finished." >>configure.log
# End of script execution
