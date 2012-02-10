#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include "vmc-typedefs.hpp"

class PositionArguments;

/**
 * Abstract base class for a lattice of any dimension.
 *
 * Any code that does not depend on the details of the lattice can use this
 * interface directly.
 */
class Lattice
{
public:
    virtual ~Lattice (void)
        {
        }

    /**
     * Returns the total number of sites on the lattice
     */
    unsigned int total_sites (void) const
        {
            return m_total_sites;
        }

    /**
     * Virtual function, called by plan_particle_move_to_nearby_empty_site()
     *
     * (plan_particle_move_to_nearby_empty_site() in random-move.hpp provides a
     * more natural API, but since the implementation depends on the details of
     * the underlying lattice, our only choice is to implement it as a virtual
     * function, which we do here.)
     *
     * @see plan_particle_move_to_nearby_empty_site()
     */
    virtual unsigned int plan_particle_move_to_nearby_empty_site_virtual (unsigned int particle, const PositionArguments &r, rng_class &rng) const = 0;

protected:
    Lattice (unsigned int total_sites)
        : m_total_sites(total_sites)
        {
        }

    const unsigned int m_total_sites;
};

#endif
