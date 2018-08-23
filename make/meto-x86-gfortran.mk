FC = gfortran
FCFLAGS_EXTRA = -ffree-line-length-none -g -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow
FCFLAGS_OPENMP = -fopenmp
PLATFORM = meto-x86-gfortran
include Makefile