#! /usr/bin/env bash

# C compiler
CC=`which icpc`; CFLAGS=" -O3 -fast -Wall"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -g -Wall"
fi
# use next line on tranditional system:
#CFLAGS="$CFLAGS `gsl-config --cflags --libs`"
# use next line on newest Ubuntu
  CFLAGS="$CFLAGS `gsl-config --cflags --libs` -B/usr/lib/i386-linux-gnu -I/usr/include/i386-linux-gnu"

# FORTRAN compiler
FC=`which ifort`; FFLAGS=" -O3 -fast -heap-arrays -cpp"
if [ "$FC" == "" ]; then
  FC=`which gfortran`; FFLAGS=" -O3 -cpp"
fi

$CC main.cpp arsenal.cpp ParameterReader.cpp RandomVariable1DArray.cpp RandomVariable2DArray.cpp RandomVariable.cpp NBD.cpp TableFunction.cpp emissionfunction.cpp Table.cpp readindata.cpp -o iSS.e $CFLAGS
