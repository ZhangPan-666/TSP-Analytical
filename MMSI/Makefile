SHELL := /bin/bash
.DEFAULT_GOAL := all
PROFILE ?= intel
ifeq ($(PROFILE),intel)
CXX = icpc
MPICXX = mpiicpc
FPFLAGS = -fp-model strict -no-ftz
OMPFLAG = -qopenmp
else ifeq ($(PROFILE),intel-llvm)
CXX = icpx
MPICXX = mpiicpx
FPFLAGS = -fp-model=precise -fno-fast-math -ffp-contract=off
OMPFLAG = -qopenmp
else ifeq ($(PROFILE),gnu)
CXX = g++
MPICXX = mpicxx
FPFLAGS = -fno-fast-math -ffp-contract=off
OMPFLAG = -fopenmp
else
$(error PROFILE must be intel, intel-llvm or gnu)
endif

CXXFLAGS ?= -O2
PROGRAMS := Work1_MassShellIntegration Work1_MagnetizationShellIntegration \
 Work1_MassShellIntegrationOpenMP Work1_MagnetizationShellIntegrationOpenMP \
 Work1_MassShellIntegrationMPI Work1_MagnetizationShellIntegrationMPI \
 Work2_MassShellIntegrationVzMPI Work2_MagnetizationShellIntegrationVectorMPI

HEADERS := _TSPMMSI.h
MPI_PROGRAMS := $(filter %MPI,$(PROGRAMS))
OMP_PROGRAMS := $(filter %OpenMP,$(PROGRAMS))
SERIAL_PROGRAMS := $(filter-out $(MPI_PROGRAMS) $(OMP_PROGRAMS),$(PROGRAMS))
DRIVER = $(CXX)
MODEFLAGS =

.PHONY: all mpi serial openmp clean
all: mpi
mpi: $(MPI_PROGRAMS)
serial: $(SERIAL_PROGRAMS)
openmp: $(OMP_PROGRAMS)
$(MPI_PROGRAMS): DRIVER = $(MPICXX)
$(OMP_PROGRAMS): MODEFLAGS = $(OMPFLAG)

$(PROGRAMS): %: %.cpp $(HEADERS) Makefile
	$(DRIVER) -I. $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(FPFLAGS) $(MODEFLAGS) "$<" -pthread -o "$@" $(LDFLAGS)

clean:
	rm -f $(PROGRAMS)
