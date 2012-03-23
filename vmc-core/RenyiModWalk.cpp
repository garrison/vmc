#include <cmath>

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

#include "RenyiModWalk.hpp"
#include "random-move.hpp"

RenyiModWalk::RenyiModWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const boost::shared_ptr<WavefunctionAmplitude> &wf_copy)
    : phialpha1(wf),
      phialpha2(wf_copy),
      transition_copy_in_progress(-1)
{
#ifndef BOOST_DISABLE_ASSERTS
    const PositionArguments &r1 = phialpha1->get_positions();
    const PositionArguments &r2 = phialpha2->get_positions();

    BOOST_ASSERT(&phialpha1->get_lattice() == &phialpha2->get_lattice());
    BOOST_ASSERT(r1.get_N_species() == r2.get_N_species());
    for (unsigned int i = 0; i < r1.get_N_species(); ++i)
        BOOST_ASSERT(r1.get_N_filled(i) == r2.get_N_filled(i));
    BOOST_ASSERT(r1.get_N_sites() == r2.get_N_sites());
    // there's no way to assert it, but we also assume they have precisely the
    // same orbitals too.  In fact, it might be useful to make a function that
    // asserts two wave functions are identical except for the particle
    // positions ...
#endif
}

probability_t RenyiModWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(transition_copy_in_progress <= 0);

    // decide which copy of the system to attempt to move
    boost::uniform_smallint<> copy_distribution(1, 2);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > copy_gen(rng, copy_distribution);
    transition_copy_in_progress = copy_gen();

    boost::shared_ptr<WavefunctionAmplitude> &phialpha = (transition_copy_in_progress == 1) ? phialpha1 : phialpha2;
    const PositionArguments &r = phialpha->get_positions();

    // plan the move
    chosen_particle = choose_random_particle(r, rng);
    const unsigned int particle_destination = plan_particle_move_to_nearby_empty_site(chosen_particle, r, phialpha->get_lattice(), rng);

#if 0
    const unsigned int N_sites = phialpha1->get_lattice().total_sites();
    std::cerr << "species " << chosen_particle.species << std::endl;

    const PositionArguments &r1 = phialpha1->get_positions();
    for (unsigned int k = 0; k < N_sites; ++k)
        std::cerr << (transition_copy_in_progress == 1 && r1[chosen_particle] == k ? '$' : (r1.is_occupied(k, chosen_particle.species) ? '*' : '-'));
    std::cerr << std::endl;

    const PositionArguments &r2 = phialpha2->get_positions();
    for (unsigned int k = 0; k < N_sites; ++k)
        std::cerr << (transition_copy_in_progress == 2 && r2[chosen_particle] == k ? '$' : (r2.is_occupied(k, chosen_particle.species) ? '*' : '-'));
    std::cerr << std::endl << std::endl;
#endif

    const amplitude_t old_phialpha_psi = phialpha->psi();

    // update phialpha using copy-on-write
    if (!phialpha.unique())
        phialpha = phialpha->clone();
    BOOST_ASSERT(phialpha.unique());

    phialpha->move_particle(chosen_particle, particle_destination);

    // return a probability
    return std::norm(phialpha->psi() / old_phialpha_psi);
}

void RenyiModWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_copy_in_progress > 0);

#if defined(DEBUG_VMC_RENYI_MOD_WALK) || defined(DEBUG_VMC_ALL)
    const PositionArguments &r1 = phialpha1->get_positions();
    for (unsigned int species = 0; species < r1.get_N_species(); ++species) {
        if (species != 0)
            std::cerr << "| ";
        for (unsigned int i = 0; i < r1.get_N_filled(species); ++i)
            std::cerr << r1[Particle(i, species)] << ' ';
    }
    std::cerr << std::endl;

    const PositionArguments &r2 = phialpha2->get_positions();
    for (unsigned int species = 0; species < r2.get_N_species(); ++species) {
        if (species != 0)
            std::cerr << "| ";
        for (unsigned int i = 0; i < r2.get_N_filled(species); ++i)
            std::cerr << r2[Particle(i, species)] << ' ';
    }
#endif

    ((transition_copy_in_progress == 1) ? phialpha1 : phialpha2)->finish_particle_moved_update();

    // remember what we just did, so each RenyiModMeasurement can update its
    // SwappedSystem
    transition_copy_just_completed = transition_copy_in_progress;

    transition_copy_in_progress = 0;
}
