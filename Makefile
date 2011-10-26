EIGEN3_CFLAGS = -I/usr/include/eigen3
BOOST_CFLAGS = 

COMPILE = g++ -Wall -W -O2 $(EIGEN3_CFLAGS) $(BOOST_CFLAGS) -DCAREFUL

HEADER_FILES = \
               CeperlyMatrix.hpp \
               MetropolisSimulation.hpp \
               Chain1d.hpp \
               SwappedSystem.hpp \
               Subsystem.hpp \
               random-combination.hpp \
               vmc-typedefs.hpp

all:	vmc

vmc:	main.o Chain1d.o random-combination.o
	$(COMPILE) -o vmc main.o Chain1d.o random-combination.o

main.o:	$(HEADER_FILES) main.cpp
	$(COMPILE) -c main.cpp

SimpleWaveFunction.o:	$(HEADER_FILES) SimpleWaveFunction.cpp
	$(COMPILE) -c SimpleWaveFunction.cpp

Chain1d.o:	$(HEADER_FILES) Chain1d.cpp
	$(COMPILE) -c Chain1d.cpp

random-combination.o:	$(HEADER_FILES) random-combination.cpp
	$(COMPILE) -c random-combination.cpp

clean:
	rm -f *.o vmc

.PHONY:	all clean
