#ifndef _RANDOM_FILLER_HPP
#define _RANDOM_FILLER_HPP

#include <vector>

#include "vmc-typedefs.hpp"

/**
 * Abstract base class used for random fillings of arbitrary lattices
 *
 * This class exists because some wavefunctions need to be in control of how
 * the different species of particles are filled due to a Gutzwiller
 * projection, but yet the WavefunctionAmplitude object is unaware of the
 * dimension of the lattice, and thus unaware of concepts like "rungs" and
 * "legs" (and therefore unable to attempt fillings with one particle per rung,
 * etc.)  So we can pass to the wavefunction a RandomFiller object, which is
 * instantiated as an NDRandomFiller object (see random-filling.hpp), and which
 * knows how to do a random filling in a smart way for a lattice of a given
 * dimension.
 */
class RandomFiller
{
public:
    virtual std::vector<unsigned int> some_random_filling (unsigned int N_filled, rng_class &rng) const = 0;

    virtual ~RandomFiller (void)
        {
        }
};

#endif
