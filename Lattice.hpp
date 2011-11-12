#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include <boost/assert.hpp>

#include "vmc-typedefs.hpp"

class PositionArguments;

class Lattice
// abstract base class.  subclasses should also implement a few implied
// functions; see e.g. HypercubicLattice.hpp
{
public:
    virtual ~Lattice (void)
	{
	}

    unsigned int total_sites (void) const
	{
	    return m_total_sites;
	}

    virtual unsigned int plan_particle_move_to_nearby_empty_site_virtual (unsigned int particle, const PositionArguments &r, rng_class &rng) const = 0;

protected:
    Lattice (unsigned int total_sites)
	: m_total_sites(total_sites)
	{
	}

    const unsigned int m_total_sites;

private:
    // disable default constructor
    Lattice (void);
};

#endif
