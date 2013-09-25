#ifndef _VMC_RANDOM_COMBINATION_HPP
#define _VMC_RANDOM_COMBINATION_HPP

#include <vector>

#include "vmc-typedefs.hpp"

class RandomNumberGenerator;

/**
 * Samples `r` indices from `0` to `n - 1`
 *
 * @param v the output will be placed in this vector
 * @param r desired length of `v`
 * @param n the upper bound (non-inclusive) of the returned indices
 * @param rng `RandomNumberGenerator` to use
 * @param keep if this is nonzero, it will leave the first `keep` elements of
 *        `v` in place.  These elements will not be repeated in the returned
 *        combination.
 */
extern void random_combination (std::vector<unsigned int> &v, unsigned int r, unsigned int n, RandomNumberGenerator &rng, unsigned int keep=0);

#endif
