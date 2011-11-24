#ifndef _RANDOM_MOVE_HPP
#define _RANDOM_MOVE_HPP

#include "PositionArguments.hpp"
#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

extern unsigned int choose_random_particle (const PositionArguments &r, rng_class &rng, unsigned int n_fixed_particles=0);

extern unsigned int choose_random_empty_site (const PositionArguments &r, rng_class &rng);

static inline unsigned int plan_particle_move_to_nearby_empty_site (unsigned int particle, const PositionArguments &r, const Lattice &lattice, rng_class &rng)
{
    // we want this to execute slightly different code based on the number of
    // dimensions in the lattice.  unfortunately, the only way to handle this
    // is through a virtual function in a Lattice subclass, so we call that
    // here
    return lattice.plan_particle_move_to_nearby_empty_site_virtual(particle, r, rng);
}

#endif
