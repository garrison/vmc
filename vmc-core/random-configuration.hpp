#ifndef _VMC_RANDOM_CONFIGURATION_HPP
#define _VMC_RANDOM_CONFIGURATION_HPP

#include <vector>

#include "Lattice.hpp"
#include "vmc-typedefs.hpp"
#include "Wavefunction.hpp"
#include "random-combination.hpp"

class RandomNumberGenerator;

/**
 * Constructs a "random" configuration of the lattice at a given filling.
 *
 * This function will at times perform some heuristics to try to get a nonzero
 * wavefunction.  For example, one such heuristic is to place each particle on
 * a different leg of a ladder until all are full.
 */
extern std::vector<unsigned int> some_random_configuration (unsigned int N_filled, const Lattice &lattice, RandomNumberGenerator &rng);

#endif
