#include <cmath>

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

#include "RenyiModWalk.hpp"
#include "random-move.hpp"

RenyiModWalk::RenyiModWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const std::vector<boost::shared_ptr<const Subsystem> > &subsystems, rng_class &rng)
    : phialpha1(wf),
      phialpha2(wf),
      transition_copy_in_progress(0)
{
    // FIXME: phialpha2 should be different than phialpha1; just take an
    // argument or rearrange the positions randomly
    (void) rng; // FIXME!

    swapped_system.reserve(subsystems.size());
    for (unsigned int i = 0; i < subsystems.size(); ++i) {
	swapped_system.push_back(boost::make_shared<SwappedSystem>(subsystems[i]));
	// fixme: don't do this until we've reached equilibrium
	swapped_system[i]->initialize(*phialpha1, *phialpha2);
    }
}

probability_t RenyiModWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(!transition_copy_in_progress);

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

    const PositionArguments &r1 = phialpha1->get_positions();
    for (unsigned int k = 0; k < N_sites; ++k)
	std::cerr << (transition_copy_in_progress == 1 && r1[chosen_particle] == k ? '$' : (r1.is_occupied(k) ? '*' : '-'));
    std::cerr << std::endl;

    const PositionArguments &r2 = phialpha2->get_positions();
    for (unsigned int k = 0; k < N_sites; ++k)
	std::cerr << (transition_copy_in_progress == 2 && r2[chosen_particle] == k ? '$' : (r2.is_occupied(k) ? '*' : '-'));
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
    BOOST_ASSERT(transition_copy_in_progress);

#ifdef DEBUG
    const PositionArguments &r1 = phialpha1->get_positions();
    for (unsigned int i = 0; i < r1.size(); ++i)
	std::cerr << r1[i] << ' ';
    std::cerr << std::endl;

    const PositionArguments &r2 = phialpha2->get_positions();
    for (unsigned int i = 0; i < r2.size(); ++i)
	std::cerr << r2[i] << ' ';
    std::cerr << std::endl << std::endl;
#endif

    ((transition_copy_in_progress == 1) ? phialpha1 : phialpha2)->finish_particle_moved_update();

    const int arg1 = (transition_copy_in_progress == 1) ? (int) chosen_particle : -1;
    const int arg2 = (transition_copy_in_progress == 2) ? (int) chosen_particle : -1;

    for (unsigned int i = 0; i < swapped_system.size(); ++i) {
	// copy on write
	if (!swapped_system[i].unique())
	    swapped_system[i] = boost::make_shared<SwappedSystem>(*swapped_system[i]);

	swapped_system[i]->update(arg1, arg2, *phialpha1, *phialpha2);
	swapped_system[i]->finish_update(*phialpha1, *phialpha2);
    }

    transition_copy_in_progress = 0;
}
