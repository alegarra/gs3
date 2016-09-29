#!/bin/bash
#this is a compilation file for gs_cage
# but this is NOT a makefile
# it compiles everything in a correct order
rm *.mod
rm *.o
rm gs3 
#ifort -c -g -traceback -check all  ALliball.f90
ifort -c -O3 ALliball.f90
ifort -c  -O3 *.f
ifort -c  -O3 sparse.f90
ifort -c  -O3 fspak90.f90
ifort -c pcg.f90
ifort -c -O3 -heap-arrays gs3.f90 
ifort *.o -O3 -heap-arrays -static  -o gs3
#ifort -c -g -traceback  -check all  -heap-arrays gs3.f90 
#ifort *.o -g -traceback -check all   -heap-arrays -static  -o gs3
