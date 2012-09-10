#ifndef _RANDOM_CONFIGURATION_HPP
#define _RANDOM_CONFIGURATION_HPP

#include <vector>

#include "Lattice.hpp"
#include "vmc-typedefs.hpp"
#include "Wavefunction.hpp"
#include "random-combination.hpp"

class RandomNumberGenerator;

extern std::vector<unsigned int> some_random_configuration (unsigned int N_filled, const Lattice &lattice, RandomNumberGenerator &rng);

extern bool search_for_configuration_with_nonzero_amplitude (Wavefunction::Amplitude &wf, RandomNumberGenerator &rng);

#endif
