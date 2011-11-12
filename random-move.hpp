#ifndef _RANDOM_MOVE_HPP
#define _RANDOM_MOVE_HPP

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/assert.hpp>

#include "PositionArguments.hpp"
#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

extern unsigned int choose_random_particle (const PositionArguments &r, rng_class &rng, unsigned int n_fixed_particles=0);

extern unsigned int choose_random_empty_site (const PositionArguments &r, rng_class &rng);

static inline unsigned int plan_particle_move_to_nearby_empty_site (unsigned int particle, const PositionArguments &r, const Lattice &lattice, rng_class &rng)
{
    // we want this to execute slightly different code based on the type of
    // lattice.  unfortunately, the only way to handle this is through a
    // virtual function in the Lattice subclass, so we call that here (and the
    // virtual function does nothing but call the template below)
    return lattice.plan_particle_move_to_nearby_empty_site_virtual(particle, r, rng);
}

template<class Lattice_T>
static inline unsigned int plan_particle_move_to_nearby_empty_site_template (unsigned int particle, const PositionArguments &r, const Lattice_T &lattice, rng_class &rng)
{
    unsigned int move_axis;
    if (Lattice_T::move_axes == 1) {
	move_axis = 0;
    } else {
	boost::uniform_smallint<> axis_distribution(0, Lattice_T::move_axes - 1);
	boost::variate_generator<rng_class&, boost::uniform_smallint<> > axis_gen(rng, axis_distribution);
	move_axis = axis_gen();
    }

    boost::uniform_smallint<> direction_distribution(0, 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > direction_gen(rng, direction_distribution);
    int step_direction = direction_gen() * 2 - 1;

    typename Lattice_T::Site site = lattice.site_from_index(r[particle]);
    unsigned int site_index;
    do {
	lattice.move_site(site, move_axis, step_direction);
	site_index = lattice.site_to_index(site);
    } while (r.is_occupied(site_index) && site_index != r[particle]);

    return site_index;
}

#endif
