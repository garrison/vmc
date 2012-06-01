#include <vector>

#include "WavefunctionAmplitude.hpp"
#include "random-filling.hpp"

void WavefunctionAmplitude::reset_with_random_positions (rng_class &rng)
{
    std::vector<std::vector<unsigned int> > vv;
    for (unsigned int i = 0; i < r.get_N_species(); ++i)
        vv.push_back(some_random_filling(r.get_N_filled(i), *lattice, rng));
    reset(PositionArguments(vv, lattice->total_sites()));
}
