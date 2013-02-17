#ifndef _VMC_RANDOM_COMBINATION_HPP
#define _VMC_RANDOM_COMBINATION_HPP

#include <vector>

#include "vmc-typedefs.hpp"

class RandomNumberGenerator;

extern void random_combination (std::vector<unsigned int> &v, unsigned int r, unsigned int n, RandomNumberGenerator &rng, unsigned int keep=0);

#endif
