#!/bin/bash
# this is NOT a Makefile, but it compiles everything in a correct order

if (( $# != 1 )); then
  echo "usage: $0 compiler"
  echo "with compiler being ifort or gfortran or x86_64-w64-mingw32-gfortran or ..."
  exit 1
fi
COMPILER=$1
# COMPILER=ifort
# COMPILER=gfortran
# COMPILER=x86_64-w64-mingw32-gfortran

hash $COMPILER 2>/dev/null || { echo >&2 "the required compiler '$COMPILER' is not installed"; exit 1; }

rm -f *.mod
rm -f *.o
rm -f gs3 gs3.exe

# $COMPILER -c -g -traceback -check all  ALliball.f90
$COMPILER -c -O3 ALliball.f90
$COMPILER -c -O3 *.f
$COMPILER -c -O3 sparse.f90
$COMPILER -c -O3 fspak90.f90
$COMPILER -c pcg.f90

if [ "$COMPILER" == "ifort" ]; then
  $COMPILER -c -O3 -heap-arrays gs3.f90
  # $COMPILER -c -g -traceback  -check all  -heap-arrays gs3.f90 
else
  $COMPILER -c -O3 gs3.f90
fi

if [ "$COMPILER" == "ifort" ]; then
  $COMPILER *.o -O3 -heap-arrays -static  -o gs3
  # $COMPILER *.o -g -traceback -check all   -heap-arrays -static  -o gs3
else
  $COMPILER *.o -O3 -static  -o gs3
fi

echo "the executable 'gs3' is available"
echo "if it was compiled for Windows, rename it as 'gs3.exe'"
