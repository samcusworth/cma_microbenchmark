#!-----------------------------------------------------------------------------
#! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
#! For further details please refer to the file LICENCE.original which you
#! should have received as part of this distribution.
#!-----------------------------------------------------------------------------
FC = ifort
FCFLAGS_EXTRA =  -qoverride-limits -g -traceback -check all -fpe0 -stand f08 -O2 -fp-model strict -warn all
FCFLAGS_OPENMP = -qopenmp
PLATFORM = meto-x86-intel
include Makefile