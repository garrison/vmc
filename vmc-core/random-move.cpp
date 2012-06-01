#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "random-move.hpp"

Particle choose_random_particle (const PositionArguments &r, rng_class &rng)
{
    // each particle has an equal probability of being chosen, regardless of
    // which species it is
    boost::uniform_smallint<> integer_distribution(0, r.get_N_filled_total() - 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > particle_gen(rng, integer_distribution);

    // this will be less than the total number of all particles in the system,
    // regardless of species
    unsigned int particle = particle_gen();

    // now we need to convert it into a species and particle number, which we
    // return
    for (unsigned int species = 0; ; ++species) {
        BOOST_ASSERT(species < r.get_N_species());
        if (particle < r.get_N_filled(species))
            return Particle(particle, species);
        particle -= r.get_N_filled(species);
    }
}

unsigned int choose_random_empty_site (const PositionArguments &r, unsigned int species, rng_class &rng)
{
    const unsigned int empty_sites = r.get_N_sites() - r.get_N_filled(species);
    BOOST_ASSERT(empty_sites > 0);

    boost::uniform_smallint<> empty_site_distribution(0, empty_sites - 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > empty_site_gen(rng, empty_site_distribution);

    int empty_site = empty_site_gen();
    for (unsigned int current_site = 0; ; ++current_site) {
        BOOST_ASSERT(current_site < r.get_N_sites());
        if (!r.is_occupied(current_site, species)) {
            if (empty_site-- == 0) {
                return current_site;
            }
        }
    }
}

unsigned int plan_particle_move_to_nearby_empty_site (Particle particle, const PositionArguments &r, const Lattice &lattice, rng_class &rng)
{
    boost::uniform_smallint<> method_distribution(0, 9);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > method_gen(rng, method_distribution);

    // occasionally, in order to make sure things are ergodic, we should
    // attempt to move to a random empty site in the lattice
    if (method_gen() == 0)
        return choose_random_empty_site(r, particle.species, rng);

    // otherwise, we find a "nearby" site
    BOOST_ASSERT(r.particle_is_valid(particle));

    unsigned int move_axis;
    if (lattice.move_axes_count() == 1) {
        move_axis = 0;
    } else {
        boost::uniform_smallint<> axis_distribution(0, lattice.move_axes_count() - 1);
        boost::variate_generator<rng_class&, boost::uniform_smallint<> > axis_gen(rng, axis_distribution);
        move_axis = axis_gen();
    }

    boost::uniform_smallint<> direction_distribution(0, 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > direction_gen(rng, direction_distribution);
    int step_direction = direction_gen() * 2 - 1;

    const unsigned int original_site_index = r[particle];
    LatticeSite site(lattice.site_from_index(original_site_index));
    unsigned int site_index;
    do {
        lattice.move_site(site, move_axis, step_direction);
        site_index = lattice.site_to_index(site);
    } while (r.is_occupied(site_index, particle.species) && site_index != original_site_index);

    return site_index;
}
