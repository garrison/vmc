EIGEN3_CFLAGS = 
BOOST_CFLAGS = 

INCLUDES = $(EIGEN3_CFLAGS) $(BOOST_CFLAGS)
LIBS = 

CXX = g++
COMPILER_OPTIMIZATIONS = -O3 -march=native
COMPILER_WARNINGS = -Wall -Wextra -Wno-unused-parameter
COMPILER_DEFINES = 
COMPILE = $(CXX) $(COMPILER_OPTIMIZATIONS) $(COMPILER_WARNINGS) $(COMPILER_DEFINES)

-include Makefile-vmc.local

HEADER_FILES = \
               allowed-momentum.hpp \
               array-util.hpp \
               BoundaryCondition.hpp \
               CeperleyMatrix.hpp \
               DensityDensityMeasurement.hpp \
               FilledOrbitals.hpp \
               FreeFermionWavefunctionAmplitude.hpp \
               HypercubicLattice.hpp \
               Lattice.hpp \
               LatticeRealization.hpp \
               lowest-momenta.hpp \
               Measurement.hpp \
               MetropolisSimulation.hpp \
               NDLattice.hpp \
               NullMeasurement.hpp \
               OrbitalDefinitions.hpp \
               PositionArguments.hpp \
               random-combination.hpp \
               random-move.hpp \
               RenyiModMeasurement.hpp \
               RenyiModWalk.hpp \
               RenyiSignMeasurement.hpp \
               RenyiSignWalk.hpp \
               SimpleSubsystem.hpp \
               StandardWalk.hpp \
               Subsystem.hpp \
               SwappedSystem.hpp \
               vmc-typedefs.hpp \
               WavefunctionAmplitude.hpp

SOURCES = \
	  FreeFermionWavefunctionAmplitude.cpp \
	  main.cpp \
	  PositionArguments.cpp \
	  random-combination.cpp \
	  random-move.cpp \
	  RenyiModWalk.cpp \
	  RenyiSignWalk.cpp \
	  StandardWalk.cpp \
	  SwappedSystem.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all:	vmc

vmc:	$(OBJECTS)
	$(COMPILE) $(LIBS) -o vmc $(OBJECTS)

%.o: %.cpp $(HEADER_FILES)
	$(COMPILE) $(INCLUDES) -c -o $@ $<

clean:
	rm -f *.o vmc

.PHONY:	all clean
