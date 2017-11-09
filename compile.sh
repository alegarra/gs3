#!/bin/bash
# this is NOT a makefile, but it compiles everything in a correct order

## choose one or the other
#COMPILATOR=ifort
COMPILATOR=gfortran

rm -f *.mod
rm -f *.o
rm -f gs3 

#$COMPILATOR -c -g -traceback -check all  ALliball.f90
$COMPILATOR -c -O3 ALliball.f90
$COMPILATOR -c -O3 *.f
$COMPILATOR -c -O3 sparse.f90
$COMPILATOR -c -O3 fspak90.f90
$COMPILATOR -c pcg.f90

if [ "$COMPILATOR" == "ifort" ]; then
  $COMPILATOR -c -O3 -heap-arrays gs3.f90
  #$COMPILATOR -c -g -traceback  -check all  -heap-arrays gs3.f90 
else
  $COMPILATOR -c -O3 gs3.f90
fi

if [ "$COMPILATOR" == "ifort" ]; then
  $COMPILATOR *.o -O3 -heap-arrays -static  -o gs3
  #$COMPILATOR *.o -g -traceback -check all   -heap-arrays -static  -o gs3
else
  $COMPILATOR *.o -O3 -static  -o gs3
fi
