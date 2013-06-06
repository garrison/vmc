#include <set>
#include <utility>

#include <boost/assert.hpp>

#include "BasicOperator.hpp"

BasicOperator::BasicOperator (const std::vector<SiteHop> &hopv_, const std::shared_ptr<const Lattice> &lattice_)
    : hopv(hopv_),
      lattice(lattice_)
{
    // validate everything except N_species
    BOOST_ASSERT(is_valid(hopv, *lattice, -1));
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
    const unsigned int n_dimensions = lattice.n_dimensions();

    // the unsigned int represents species
    std::set<std::pair<LatticeSite, unsigned int> > unique_pair_set;

    for (std::vector<SiteHop>::const_iterator i = hopv.begin(); i != hopv.end(); ++i) {
        // check number of dimensions for each given LatticeSite, and check the
        // species index given
        if (!(i->get_source().n_dimensions() == n_dimensions
              && i->get_destination().n_dimensions() == n_dimensions
              && i->get_species() < N_species))
            return false;

        // check that LatticeSite's are all unique for a species, up to modulo
        // of the lattice (except the number counting operator, which must be
        // unique from all others)
        if (!unique_pair_set.insert(std::make_pair(inline_enforce_boundary(i->get_source(), lattice), i->get_species())).second)
            return false;
        if (i->get_source() != i->get_destination() && !unique_pair_set.insert(std::make_pair(inline_enforce_boundary(i->get_destination(), lattice), i->get_species())).second)
            return false;
    }

    return true;
}
