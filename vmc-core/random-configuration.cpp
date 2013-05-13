#include <set>
#include <algorithm>

#include <boost/assert.hpp>

#include "random-configuration.hpp"
#include "RandomNumberGenerator.hpp"

std::vector<unsigned int> some_random_configuration (unsigned int N_filled, const Lattice &lattice, RandomNumberGenerator &rng)
{
    BOOST_ASSERT(N_filled <= lattice.total_sites());

    const unsigned int n_dimensions = lattice.n_dimensions();

    // If in more than one dimension, occasionally we want to try placing
    // particles such that they are distributed as well as possible over a
    // certain dimension (rungs, legs, etc).  For determinantal wavefunctions
    // with certain orbital configurations, this can often be necessary to find
    // a non-zero amplitude in a reasonable amount of time.
    //
    // This method could be further optimized, but it is fast enough for now
    // (and is unlikely to be the bottleneck anyway).
    if (n_dimensions > 1 && rng.random_small_uint(3) == 0) {
        std::vector<unsigned int> v;
        std::set<unsigned int> vs;
        unsigned int spread_dimension = rng.random_small_uint(n_dimensions);
        unsigned int remaining = N_filled;
        while (remaining != 0) {
            std::vector<unsigned int> spread_coordinate;
            random_combination(spread_coordinate,
                               std::min(remaining, (unsigned int) lattice.dimensions[spread_dimension]),
                               lattice.dimensions[spread_dimension], rng);
            for (unsigned int i = 0; i < spread_coordinate.size(); ++i) {
                unsigned int proposed_site_index;
                do {
                    LatticeSite proposed_site(lattice.n_dimensions());
                    proposed_site[spread_dimension] = spread_coordinate[i];
                    for (unsigned int j = 0; j < n_dimensions; ++j) {
                        if (j != spread_dimension)
                            proposed_site[j] = rng.random_small_uint(lattice.dimensions[j]);
                    }
                    proposed_site.basis_index = rng.random_small_uint(lattice.basis_indices);
                    proposed_site_index = lattice.index(proposed_site);
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
