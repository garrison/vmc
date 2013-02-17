#ifndef _VMC_BASIC_OPERATOR_HPP
#define _VMC_BASIC_OPERATOR_HPP

#include <vector>

#include <boost/shared_ptr.hpp>

#include "Lattice.hpp"

class SiteHop
{
public:
    SiteHop (const LatticeSite &source_, const LatticeSite &destination_, unsigned int species_)
        : source(source_),
          destination(destination_),
          species(species_)
        {
        }

    const LatticeSite & get_source (void) const
        {
            return source;
        }

    const LatticeSite & get_destination (void) const
        {
            return destination;
        }

    unsigned int get_species (void) const
        {
            return species;
        }

private:
    // these should be constant after initialization
    LatticeSite source, destination;
    unsigned int species;
};

/**
 * An operator that conserves particle number
 *
 * The operator must conserve particle number
 *
 * and be written like creation * annihilation * creation * annihilation (...)
 *
 * (i.e. a series of hops)
 *
 * A given site can only be mentioned once.  The only exception is to mark the
 * same site as source and destination, in which case it evaluates to 1 or 0
 * based on whether a particle is present on the site (just like the ordinary
 * number counting operator).
 *
 * Be careful putting a large operator on a small lattice, as it might wrap
 * over the boundary and no longer make sense in the way intended (or might not
 * even be valid, for that matter)
 */
class BasicOperator
{
public:
    BasicOperator (const std::vector<SiteHop> &hopv_, const boost::shared_ptr<const Lattice> &lattice_);

    static bool is_valid (const std::vector<SiteHop> &hopv, const Lattice &lattice, unsigned int N_species);

    const std::vector<SiteHop> hopv;
    const boost::shared_ptr<const Lattice> lattice;
};

#endif
