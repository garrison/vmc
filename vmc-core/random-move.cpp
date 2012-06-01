#include <boost/assert.hpp>

#include "RandomNumberGenerator.hpp"
#include "random-move.hpp"

Particle choose_random_particle (const PositionArguments &r, RandomNumberGenerator &rng)
{
    // Each particle has an equal probability of being chosen, regardless of
    // which species it is.  The following quantity will be fewer than the
    // total number of all particles in the system, regardless of species.
    unsigned int particle = rng.random_small_uint(r.get_N_filled_total());

    // now we need to convert it into a species and particle number, which we
    // return
    for (unsigned int species = 0; ; ++species) {
        BOOST_ASSERT(species < r.get_N_species());
        if (particle < r.get_N_filled(species))
            return Particle(particle, species);
        particle -= r.get_N_filled(species);
    }
}

unsigned int choose_random_empty_site (const PositionArguments &r, unsigned int species, RandomNumberGenerator &rng)
{
    const unsigned int empty_sites = r.get_N_sites() - r.get_N_filled(species);
    BOOST_ASSERT(empty_sites > 0);

    int empty_site = rng.random_small_uint(empty_sites);
    for (unsigned int current_site = 0; ; ++current_site) {
        BOOST_ASSERT(current_site < r.get_N_sites());
        if (!r.is_occupied(current_site, species)) {
            if (empty_site-- == 0) {
                return current_site;
            }
        }
    }
}

unsigned int plan_particle_move_to_nearby_empty_site (Particle particle, const PositionArguments &r, const Lattice &lattice, RandomNumberGenerator &rng)
{
    // occasionally, in order to make sure things are ergodic, we should
    // attempt to move to a random empty site in the lattice
    if (rng.random_small_uint(10) == 0)
        return choose_random_empty_site(r, particle.species, rng);

    // otherwise, we find a "nearby" site
    BOOST_ASSERT(r.particle_is_valid(particle));

    unsigned int move_axis;
    if (lattice.move_axes_count() == 1) {
        move_axis = 0;
    } else {
        move_axis = rng.random_small_uint(lattice.move_axes_count());
    }

    // will be either +1 or -1
    int step_direction = rng.random_small_uint(2) * 2 - 1;

    const unsigned int original_site_index = r[particle];
    LatticeSite site(lattice.site_from_index(original_site_index));
    unsigned int site_index;
    do {
        lattice.move_site(site, move_axis, step_direction);
        site_index = lattice.site_to_index(site);
    } while (r.is_occupied(site_index, particle.species) && site_index != original_site_index);

    return site_index;
}
