#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "random-move.hpp"

unsigned int choose_random_particle (const PositionArguments &r, rng_class &rng, unsigned int n_fixed_particles)
{
    BOOST_ASSERT(n_fixed_particles < r.get_N_filled());
    boost::uniform_smallint<> integer_distribution(n_fixed_particles, r.get_N_filled() - 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > particle_gen(rng, integer_distribution);

    return particle_gen();
}

unsigned int choose_random_empty_site (const PositionArguments &r, rng_class &rng)
{
    unsigned int empty_sites = r.get_N_sites() - r.get_N_filled();
    BOOST_ASSERT(empty_sites > 0);

    boost::uniform_smallint<> empty_site_distribution(0, empty_sites - 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > empty_site_gen(rng, empty_site_distribution);

    int empty_site = empty_site_gen();
    for (unsigned int current_site = 0; ; ++current_site) {
        BOOST_ASSERT(current_site < r.get_N_sites());
        if (!r.is_occupied(current_site)) {
            if (empty_site-- == 0) {
                return current_site;
            }
        }
    }
}
