#ifndef _VMC_RANDOM_CONFIGURATION_HPP
#define _VMC_RANDOM_CONFIGURATION_HPP

#include <vector>

#include "Lattice.hpp"
#include "vmc-typedefs.hpp"
#include "Wavefunction.hpp"
#include "random-combination.hpp"

class RandomNumberGenerator;

extern std::vector<unsigned int> some_random_configuration (unsigned int N_filled, const Lattice &lattice, RandomNumberGenerator &rng);

#endif
