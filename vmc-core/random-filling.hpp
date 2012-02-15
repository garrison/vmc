#ifndef _RANDOM_FILLING_HPP
#define _RANDOM_FILLING_HPP

#include "NDLattice.hpp"
#include "vmc-typedefs.hpp"
#include "random-combination.hpp"
#include "WavefunctionAmplitude.hpp"

// fixme: there is currently no reason for this to be a template
template <unsigned int DIM>
PositionArguments some_random_filling (unsigned int N_filled, const NDLattice<DIM> &lattice, rng_class &rng)
{
    std::vector<unsigned int> v;
    random_combination(v, N_filled, lattice.total_sites(), rng);
    return PositionArguments(v, lattice.total_sites());
}

// fixme: there is currently no reason for this to be a template
template <unsigned int DIM>
bool search_for_filling_with_nonzero_amplitude (WavefunctionAmplitude &wf, const NDLattice<DIM> &lattice, rng_class &rng)
{
    unsigned int attempts = 1; // assume that one attempt has already been completed
    while (wf.psi() == amplitude_t(0)) {
        if (attempts++ == 10000)
            return false;
        wf.reset(some_random_filling<DIM>(wf.get_positions().get_N_filled(), lattice, rng));
    }
    return true;
}

#endif
