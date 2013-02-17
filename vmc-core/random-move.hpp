#ifndef _VMC_RANDOM_MOVE_HPP
#define _VMC_RANDOM_MOVE_HPP

#include "PositionArguments.hpp"
#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

class RandomNumberGenerator;

/**
 * Chooses a random particle index.
 *
 * @param r PositionArguments representing the particles in the system
 * @param rng random number generator
 * @return the chosen particle
 */
extern Particle choose_random_particle (const PositionArguments &r, RandomNumberGenerator &rng);

/**
 * Chooses a random empty site
 *
 * NOTE: "empty" here means that there is not already a particle of the given
 * species on the site.  An "empty" site may in fact contain particle(s) of a
 * different species.
 *
 * @param r PositionArguments representing the particles in the system
 * @param rng random number generator
 * @return index of chosen empty site, which is a number less than
 *         r.get_N_sites()
 */
extern unsigned int choose_random_empty_site (const PositionArguments &r, unsigned int species, RandomNumberGenerator &rng);

/**
 * Returns an empty site index that is "near" a given particle.
 *
 * The method used to return a nearby site should satisfy balance, so it can be
 * used to plan moves in Monte Carlo simulations.
 *
 * On occasion, this function will return a site that is far away, just to make
 * sure things are ergodic in the simulation.
 *
 * NOTE: "empty" here means that there is not already a particle of the given
 * species on the site.  An "empty" site may in fact contain particle(s) of a
 * different species.
 *
 * @param particle index of particle to move
 * @param r PositionArguments representing the particles in the system
 * @param lattice the lattice itself
 * @param rng random number generator
 * @return index of chosen nearby empty site, which is a number less than
 *         r.get_N_sites()
 */
extern unsigned int plan_particle_move_to_nearby_empty_site (Particle particle, const PositionArguments &r, const Lattice &lattice, RandomNumberGenerator &rng);

#endif
