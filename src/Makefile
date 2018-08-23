ifndef PLATFORM
$(error Platform file not loaded, re-run as "make -f" providing a file from \
the "../make" subdirectory as an argument)
endif

FCFLAGS = ${FCFLAGS_EXTRA}

BUILD_OPENMP ?= true
ifeq (${BUILD_OPENMP}, true)
FCFLAGS += ${FCFLAGS_OPENMP}
endif

SRCS = $(patsubst %.F90, %.o, $(wildcard *.F90))

PROGRAM = prog.exe

$(PROGRAM): cma_test_algorithm_mod_psy.f90 $(SRCS)
	${FC} ${FCFLAGS} -o $@ $^

# this module has an explicit dependencies
columnwise_op_mul_kernel_mod.o: argument_mod.o constants_mod.o kernel_mod.o

%.o: %.F90
	${FC} ${FCFLAGS} -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.o *.mod