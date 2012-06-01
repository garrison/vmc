#ifndef _RANDOM_FILLING_HPP
#define _RANDOM_FILLING_HPP

#include <vector>

#include "Lattice.hpp"
#include "vmc-typedefs.hpp"
#include "WavefunctionAmplitude.hpp"
#include "random-combination.hpp"

extern std::vector<unsigned int> some_random_filling (unsigned int N_filled, const Lattice &lattice, rng_class &rng);

extern bool search_for_filling_with_nonzero_amplitude (WavefunctionAmplitude &wf, rng_class &rng);

#endif
