#include <set>
#include <utility>
#include <algorithm>
#include <cassert>

#include "BasicOperator.hpp"

BasicOperator::BasicOperator (const std::vector<SiteHop> &hopv_, const std::shared_ptr<const Lattice> &lattice_)
    : hopv(hopv_),
      lattice(lattice_)
{
    // validate everything except N_species
    assert(is_valid(hopv, *lattice, -1));
}

static inline LatticeSite inline_enforce_boundary (const LatticeSite &site, const Lattice &lattice)
{
    LatticeSite rv(site);
    // ignore the phase change (which is not necessary/meaningful here)
    lattice.enforce_periodic_boundary(rv);
    return rv;
}

bool BasicOperator::is_valid (const std::vector<SiteHop> &hopv, const Lattice &lattice, unsigned int N_species)
{
    return (min_required_species(hopv, lattice) <= N_species);
}

unsigned int BasicOperator::min_required_species (const std::vector<SiteHop> &hopv, const Lattice &lattice)
{
    const unsigned int n_dimensions = lattice.n_dimensions();

    if (hopv.empty())
        return 0; // do not allow empty operators

    unsigned int max_species_index = 0;

    // the unsigned int represents species
    std::set<std::pair<LatticeSite, unsigned int> > unique_pair_set;

    for (std::vector<SiteHop>::const_iterator i = hopv.begin(); i != hopv.end(); ++i) {
        // check number of dimensions for each given LatticeSite
        if (i->get_source().n_dimensions() != n_dimensions
            || i->get_destination().n_dimensions() != n_dimensions)
            return 0;

        // keep track of the highest species index referenced by the sitehops
        max_species_index = std::max(max_species_index, i->get_species());

        // check that LatticeSite's are all unique for a species, up to modulo
        // of the lattice (except the number counting operator, which must be
        // unique from all others)
        if (!unique_pair_set.insert(std::make_pair(inline_enforce_boundary(i->get_source(), lattice), i->get_species())).second)
            return 0;
        if (i->get_source() != i->get_destination() && !unique_pair_set.insert(std::make_pair(inline_enforce_boundary(i->get_destination(), lattice), i->get_species())).second)
            return 0;
    }

    return max_species_index + 1;
}
