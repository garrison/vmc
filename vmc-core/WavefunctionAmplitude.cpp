#include <vector>

#include "WavefunctionAmplitude.hpp"

void WavefunctionAmplitude::reset_with_filler (const RandomFiller &filler, rng_class &rng)
{
    std::vector<std::vector<unsigned int> > vv;
    for (unsigned int i = 0; i < r.get_N_species(); ++i)
        vv.push_back(filler.some_random_filling(r.get_N_filled(i), rng));
    reset(PositionArguments(vv, lattice->total_sites()));
}
