EIGEN3_CFLAGS = 
BOOST_CFLAGS = 
JSONCPP_CFLAGS = 
JSONCPP_LIBS = 

INCLUDES = $(EIGEN3_CFLAGS) $(BOOST_CFLAGS) $(JSONCPP_CFLAGS) $(EXTRA_CFLAGS)
LIBS = $(JSONCPP_LIBS) $(EXTRA_LIBS)

CXX = g++
COMPILER_OPTIMIZATIONS = -O3 #-march=native
COMPILER_WARNINGS = -Wall #-Wextra -Wno-unused-parameter
COMPILER_DEFINES = #-DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG
COMPILE = $(CXX) $(COMPILER_OPTIMIZATIONS) $(COMPILER_WARNINGS) $(COMPILER_DEFINES)

DOXYGEN = doxygen

-include Makefile-vmc.local

HEADER_FILES = \
               allowed-momentum.hpp \
               BoundaryCondition.hpp \
               CeperleyMatrix.hpp \
               DBLWavefunctionAmplitude.hpp \
               DMetalWavefunctionAmplitude.hpp \
               DensityDensityMeasurement.hpp \
               FilledOrbitals.hpp \
               FreeFermionWavefunctionAmplitude.hpp \
               GreenMeasurement.hpp \
               HypercubicLattice.hpp \
               Lattice.hpp \
               LatticeRealization.hpp \
               lowest-momenta.hpp \
               Measurement.hpp \
               MetropolisSimulation.hpp \
               NDLattice.hpp \
               OrbitalDefinitions.hpp \
               PositionArguments.hpp \
               random-combination.hpp \
               random-move.hpp \
               RenyiModMeasurement.hpp \
               RenyiModWalk.hpp \
               RenyiSignMeasurement.hpp \
               RenyiSignWalk.hpp \
               RunInformation.hpp \
               SimpleSubsystem.hpp \
               StandardWalk.hpp \
               Subsystem.hpp \
               SwappedSystem.hpp \
               vmc-math-utils.hpp \
               vmc-typedefs.hpp \
               WavefunctionAmplitude.hpp

SOURCES = \
          CeperleyMatrix.cpp \
          DBLWavefunctionAmplitude.cpp \
          DMetalWavefunctionAmplitude.cpp \
          FreeFermionWavefunctionAmplitude.cpp \
          main.cpp \
          PositionArguments.cpp \
          random-combination.cpp \
          random-move.cpp \
          RenyiModWalk.cpp \
          RenyiSignWalk.cpp \
          RunInformation.cpp \
          StandardWalk.cpp \
          SwappedSystem.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all:	vmc-core

vmc-core: $(OBJECTS)
	$(COMPILE) $(LIBS) -o vmc-core $(OBJECTS)

%.o: %.cpp $(HEADER_FILES) Makefile Makefile-vmc.local
	$(COMPILE) $(INCLUDES) -DVMC_COMPILER_COMMANDLINE="$(COMPILE)" -c -o $@ $<

docs: Doxyfile $(HEADER_FILES) $(SOURCES) clean_docs
	mkdir -p docs/
	$(DOXYGEN)

clean_docs:
	rm -rf docs/generated/

clean: clean_docs
	rm -f *.o vmc

.PHONY:	all clean docs clean_docs
