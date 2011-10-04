EIGEN3_INCLUDE_PATH = /usr/include/eigen3

COMPILE = g++ -Wall -W -O2 -DDEBUG -I$(EIGEN3_INCLUDE_PATH)

HEADER_FILES = \
               CeperlyMatrix.hpp \
               MetropolisSimulation.hpp \
               Chain1d.hpp \
               vmc-typedefs.hpp

all:	vmc

vmc:	main.o Chain1d.o
	$(COMPILE) -o vmc main.o Chain1d.o

main.o:	$(HEADER_FILES) main.cpp
	$(COMPILE) -c main.cpp

SimpleWaveFunction.o:	$(HEADER_FILES) SimpleWaveFunction.cpp
	$(COMPILE) -c SimpleWaveFunction.cpp

Chain1d.o:	$(HEADER_FILES) Chain1d.cpp
	$(COMPILE) -c Chain1d.cpp

clean:
	rm -f *.o vmc

.PHONY:	all clean
