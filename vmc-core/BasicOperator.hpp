#ifndef _VMC_BASIC_OPERATOR_HPP
#define _VMC_BASIC_OPERATOR_HPP

#include <vector>
#include <memory>

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
 * A "basic" operator
 *
 * The operator must conserve the particle number of each species and be
 * written as a product of (creation * annihilation) factors (i.e. a series of
 * hops, each of which acts on a single species).
 *
 * A given site + species combination can only be mentioned once in the string.
 * The only exception is to mark the same site as source and destination, in
 * which case it evaluates to 1 or 0 based on whether a particle is present on
 * the site (just like the ordinary number counting operator).
 *
 * Be careful putting a large operator on a small lattice, as it might wrap
 * over the boundary and no longer make sense in the way intended (or might not
 * even be valid, for that matter)
 */
class BasicOperator
{
public:
    BasicOperator (const std::vector<SiteHop> &hopv_, const std::shared_ptr<const Lattice> &lattice_);

    /**
     * Tests whether the sitehops are valid on a lattice with a given number of species
     */
    static bool is_valid (const std::vector<SiteHop> &hopv, const Lattice &lattice, unsigned int N_species);

    /**
     * Returns the number of species a wavefunction must have for the operator
     * to be valid.  If the operator is invalid on the given lattice, this
     * function returns zero.
     */
    static unsigned int min_required_species (const std::vector<SiteHop> &hopv, const Lattice &lattice);

    const std::vector<SiteHop> hopv;
    const std::shared_ptr<const Lattice> lattice;
};

#endif
