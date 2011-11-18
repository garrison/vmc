EIGEN3_CFLAGS = -I/usr/include/eigen3
BOOST_CFLAGS = 

INCLUDES = $(EIGEN3_CFLAGS) $(BOOST_CFLAGS)
LIBS = 

COMPILE = g++ -Wall -W -O2 $(INCLUDES)

HEADER_FILES = \
               BoundaryCondition.hpp \
               CeperleyMatrix.hpp \
               DensityDensityMeasurement.hpp \
               Chain1dOrbitals.hpp \
               FreeFermionWavefunctionAmplitude.hpp \
               HypercubicLattice.hpp \
               HypercubicSubsystem.hpp \
               Lattice.hpp \
               Measurement.hpp \
               MetropolisSimulation.hpp \
               NullMeasurement.hpp \
               OrbitalDefinitions.hpp \
               PositionArguments.hpp \
               random-combination.hpp \
               random-move.hpp \
               RenyiModMeasurement.hpp \
               RenyiModWalk.hpp \
               RenyiSignMeasurement.hpp \
               RenyiSignWalk.hpp \
               StandardWalk.hpp \
               Subsystem.hpp \
               SwappedSystem.hpp \
               vmc-typedefs.hpp \
               WavefunctionAmplitude.hpp

SOURCES = \
	  Chain1dOrbitals.cpp \
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
	$(COMPILE) -c -o $@ $<

clean:
	rm -f *.o vmc

.PHONY:	all clean
