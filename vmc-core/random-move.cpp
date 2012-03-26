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

    // we want this to execute slightly different code based on the number of
    // dimensions in the lattice.  unfortunately, the only way to handle this
    // is through a virtual function in a Lattice subclass, so we call that
    // here
    return lattice.plan_particle_move_to_nearby_empty_site_virtual(particle, r, rng);
}
