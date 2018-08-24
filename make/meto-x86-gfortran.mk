#!-----------------------------------------------------------------------------
#! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
#! For further details please refer to the file LICENCE.original which you
#! should have received as part of this distribution.
#!-----------------------------------------------------------------------------
FC = gfortran
FCFLAGS_EXTRA = -ffree-line-length-none -g -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow
FCFLAGS_OPENMP = -fopenmp
PLATFORM = meto-x86-gfortran
include Makefile