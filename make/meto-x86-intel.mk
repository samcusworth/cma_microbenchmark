FC = ifort
#  -qopenmp
FCFLAGS_EXTRA =  -qoverride-limits -g -traceback -check all -fpe0 -stand f08 -O2 -fp-model strict -warn all
FCFLAGS_OPENMP = -qopenmp
PLATFORM = meto-x86-intel
include Makefile