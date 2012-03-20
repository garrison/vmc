#ifndef _RANDOM_MOVE_HPP
#define _RANDOM_MOVE_HPP

#include "PositionArguments.hpp"
#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

/**
 * Chooses a random particle index.
 *
 * @param r PositionArguments representing the particles in the system
 * @param rng random number generator
 * @return the chosen particle
 */
extern Particle choose_random_particle (const PositionArguments &r, rng_class &rng);

/**
 * Chooses a random empty site
 *
 * @param r PositionArguments representing the particles in the system
 * @param rng random number generator
 * @return index of chosen empty site, which is a number less than
 *         r.get_N_sites()
 */
extern unsigned int choose_random_empty_site (const PositionArguments &r, unsigned int species, rng_class &rng);

/**
 * Returns an empty site index that is "near" a given particle.
 *
 * The method used to return a nearby site should satisfy balance, so it can be
 * used to plan moves in Monte Carlo simulations.
 *
 * @param particle index of particle to move
 * @param r PositionArguments representing the particles in the system
 * @param lattice the lattice itself
 * @param rng random number generator
 * @return index of chosen nearby empty site, which is a number less than
 *         r.get_N_sites()
 */
static inline unsigned int plan_particle_move_to_nearby_empty_site (Particle particle, const PositionArguments &r, const Lattice &lattice, rng_class &rng)
{
    // we want this to execute slightly different code based on the number of
    // dimensions in the lattice.  unfortunately, the only way to handle this
    // is through a virtual function in a Lattice subclass, so we call that
    // here
    return lattice.plan_particle_move_to_nearby_empty_site_virtual(particle, r, rng);
}

#endif
