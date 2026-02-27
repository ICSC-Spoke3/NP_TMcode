#!/bin/bash

# Script settings
l_start=2 # first order to test
l_end=20 # last order to test
wavelength="1.000E-06" # wavelength to test (use your shortest one)
dielec_real="2.6015926E+00" # real part of dielectric function
dielec_imag="3.8710688E-04" # imaginary part of dielectric function
radius=$1 # radius is read from command line
ompthreads="" # OMP threads to use (default is "1,NPROC")
outputfile="convergence.csv"

# Script execution starts here
if [ "x${ompthreads}" = "x" ]; then
    NPROC=$(nproc)
    ompthreads="1,${NPROC}" # OMP threads settings
fi
if [ "x${NPTMDIR}" = "x" ]; then
    echo "ERROR: environment variable NPTMDIR is not defined."
    echo "       Use \"export NPTMDIR=PATH_TO_NPTMCODE_INSTALL_DIR\","
    echo "       where the installation directory is the one containing"
    echo "       the COPYING license file."
    exit 1
fi
li=${l_start}
echo "#Wavelength,ScaSec,AbsSec,ExtSec,LM" > ${outputfile}
while [ $((li)) -le ${l_end} ]
do
    fold_index=$(python -c "print('{0:02d}'.format($li))")
    str_index=$(python -c "print('{0:2d}'.format($li))")
    folder="li_${fold_index}"
    echo $li
    mkdir -p $folder
    cd $folder
    cat > DEDFB <<EOF
   1  0
 1.0000000E+00 3.0000000E+15 1.0000000E+00  0    1    0   3
1.000E-06
    1
  1     ${radius}
 1.00000E+00
(${dielec_real},${dielec_imag})
0
EOF
    cat > DSPH <<EOF
    1   ${str_index}    0  149  300    0
 0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 1
0
USE_DYN_ORDERS=0
EOF
    OMP_NUM_THREADS=${ompthreads} ${NPTMDIR}/build/sphere/np_sphere DEDFB DSPH .
    ${NPTMDIR}/src/scripts/parse_output.py --app=SPH --in c_OSPH --out convergence
    cd ..
    str_line=$(tail -1 ${folder}/convergence_ics.csv)
    echo "${str_line},${fold_index}" >> ${outputfile}
    rm -rf ${folder}
    li=$((li+1))
done
