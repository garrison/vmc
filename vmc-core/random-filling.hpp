#ifndef _RANDOM_FILLING_HPP
#define _RANDOM_FILLING_HPP

#include <vector>
#include <set>
#include <algorithm>

#include <boost/assert.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include "Lattice.hpp"
#include "vmc-typedefs.hpp"
#include "random-combination.hpp"
#include "RandomFiller.hpp"
#include "WavefunctionAmplitude.hpp"

// fixme: do we really want this to be here?  if we move it, remove the boost includes above.
static inline unsigned int random_small_uint (unsigned int min, unsigned int max, rng_class &rng)
{
    boost::uniform_smallint<> distribution(min, max);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > generator(rng, distribution);
    return generator();
}

template <unsigned int DIM>
std::vector<unsigned int> some_random_filling (unsigned int N_filled, const NDLattice<DIM> &lattice, rng_class &rng)
{
    BOOST_ASSERT(N_filled <= lattice.total_sites());

    // If in more than one dimension, occasionally we want to try filling
    // particles such that they are distributed as well as possible over a
    // certain dimension (rungs, legs, etc).  For determinantal wavefunctions
    // with certain orbital configurations, this can often be necessary to find
    // a non-zero amplitude in a reasonable amount of time.
    //
    // This method could be further optimized, but it is fast enough for now
    // (and is unlikely to be the bottleneck anyway).
    if (DIM > 1 && random_small_uint(0, 2, rng) == 0) {
        std::vector<unsigned int> v;
        std::set<unsigned int> vs;
        unsigned int spread_dimension = random_small_uint(0, DIM - 1, rng);
        unsigned int remaining = N_filled;
        while (remaining != 0) {
            std::vector<unsigned int> spread_coordinate;
            random_combination(spread_coordinate,
                               std::min(remaining, (unsigned int) lattice.length[spread_dimension]),
                               lattice.length[spread_dimension], rng);
            for (unsigned int i = 0; i < spread_coordinate.size(); ++i) {
                unsigned int proposed_site_index;
                do {
                    typename NDLattice<DIM>::Site proposed_site;
                    proposed_site[spread_dimension] = spread_coordinate[i];
                    for (unsigned int j = 0; j < DIM; ++j) {
                        if (j != spread_dimension)
                            proposed_site[j] = random_small_uint(0, lattice.length[j] - 1, rng);
                    }
                    proposed_site.basis_index = random_small_uint(0, lattice.basis_indices - 1, rng);
                    proposed_site_index = lattice.site_to_index(proposed_site);
                } while (!vs.insert(proposed_site_index).second); // try again until successful
                v.push_back(proposed_site_index);
            }
            remaining -= spread_coordinate.size();
        }
        return v;
    }

    // otherwise, fall back to just blindly choosing a random combination of sites
    std::vector<unsigned int> v;
    random_combination(v, N_filled, lattice.total_sites(), rng);
    return v;
}

/**
 * N-dimensional random filler
 *
 * @see RandomFiller
 */
template <unsigned int DIM>
class NDRandomFiller : public RandomFiller
{
public:
    NDRandomFiller (const NDLattice<DIM> &lattice_)
        : lattice(lattice_)
        {
        }

    std::vector<unsigned int> some_random_filling (unsigned int N_filled, rng_class &rng) const
        {
            return ::some_random_filling<DIM>(N_filled, lattice, rng);
        }

private:
    const NDLattice<DIM> &lattice;
};

template <unsigned int DIM>
bool search_for_filling_with_nonzero_amplitude (WavefunctionAmplitude &wf, const NDLattice<DIM> &lattice, rng_class &rng)
{
    NDRandomFiller<DIM> filler(lattice);
    unsigned int attempts = 1; // assume that one attempt has already been completed
    while (wf.psi() == amplitude_t(0)) {
        if (attempts++ == 1000000)
            return false;
        wf.reset_with_filler(filler, rng);
    }
    return true;
}

#endif
