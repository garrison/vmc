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
COMPILE = $(CXX) $(COMPILER_OPTIMIZATIONS) $(COMPILER_WARNINGS)

DOXYGEN = doxygen

-include Makefile-vmc.local

HEADER_FILES = \
               BaseSwapPossibleWalk.hpp \
               BinnedEstimate.hpp \
               BoundaryCondition.hpp \
               CeperleyMatrix.hpp \
               DBLWavefunctionAmplitude.hpp \
               DMetalWavefunctionAmplitude.hpp \
               FreeFermionWavefunctionAmplitude.hpp \
               Lattice.hpp \
               lw_vector.hpp \
               Measurement.hpp \
               MetropolisSimulation.hpp \
               Move.hpp \
               OperatorMeasurement.hpp \
               OrbitalDefinitions.hpp \
               ParticleOperator.hpp \
               PositionArguments.hpp \
               RandomNumberGenerator.hpp \
               random-combination.hpp \
               random-configuration.hpp \
               random-move.hpp \
               RenyiModPossibleMeasurement.hpp \
               RenyiModPossibleWalk.hpp \
               RenyiSignMeasurement.hpp \
               RenyiSignWalk.hpp \
               RunInformation.hpp \
               RunningEstimate.hpp \
               RVBWavefunctionAmplitude.hpp \
               SimpleSubsystem.hpp \
               StandardWalk.hpp \
               Subsystem.hpp \
               SubsystemOccupationNumberProbabilityMeasurement.hpp \
               SwappedSystem.hpp \
               vmc-math-utils.hpp \
               vmc-typedefs.hpp \
               Walk.hpp \
               WavefunctionAmplitude.hpp

SOURCES = \
          BaseSwapPossibleWalk.cpp \
          CeperleyMatrix.cpp \
          DBLWavefunctionAmplitude.cpp \
          DMetalWavefunctionAmplitude.cpp \
          FreeFermionWavefunctionAmplitude.cpp \
          Lattice.cpp \
          MetropolisSimulation.cpp \
          Move.cpp \
          OperatorMeasurement.cpp \
          ParticleOperator.cpp \
          PositionArguments.cpp \
          RandomNumberGenerator.cpp \
          random-combination.cpp \
          random-configuration.cpp \
          random-move.cpp \
          RenyiModPossibleWalk.cpp \
          RenyiSignWalk.cpp \
          RunInformation.cpp \
          RVBWavefunctionAmplitude.cpp \
          SimpleSubsystem.cpp \
          StandardWalk.cpp \
          SwappedSystem.cpp \
          WavefunctionAmplitude.cpp

OBJECTS = $(SOURCES:.cpp=.o)

all:	vmc-core

vmc-core: $(OBJECTS) vmc-core.o
	$(COMPILE) $(LIBS) -o vmc-core $(OBJECTS) vmc-core.o

%.o: %.cpp $(HEADER_FILES) Makefile Makefile-vmc.local
	$(COMPILE) $(INCLUDES) $(COMPILER_DEFINES) -DVMC_COMPILER_COMMANDLINE="$(COMPILE)" -c -o $@ $<

docs: Doxyfile $(HEADER_FILES) $(SOURCES) clean_docs
	mkdir -p docs/
	$(DOXYGEN)

clean_docs:
	rm -rf docs/generated/

clean: clean_docs
	rm -f *.o vmc

.PHONY:	all clean docs clean_docs
