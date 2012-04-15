#include <cmath>

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

#include "BaseSwapPossibleWalk.hpp"
#include "random-move.hpp"

BaseSwapPossibleWalk::BaseSwapPossibleWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const boost::shared_ptr<WavefunctionAmplitude> &wf_copy, boost::shared_ptr<const Subsystem> subsystem, bool update_swapped_system_before_accepting_)
    : phialpha1(wf),
      phialpha2(wf_copy),
      swapped_system(new SwappedSystem(subsystem)),
      update_swapped_system_before_accepting(update_swapped_system_before_accepting_),
      transition_in_progress(false)
{
#ifndef BOOST_DISABLE_ASSERTS
    BOOST_ASSERT(&phialpha1->get_lattice() == &phialpha2->get_lattice());
    BOOST_ASSERT(phialpha1->get_positions().get_N_species() == phialpha2->get_positions().get_N_species());
    for (unsigned int i = 0; i < phialpha1->get_positions().get_N_species(); ++i)
        BOOST_ASSERT(phialpha1->get_positions().get_N_filled(i) == phialpha2->get_positions().get_N_filled(i));
    BOOST_ASSERT(phialpha1->get_positions().get_N_sites() == phialpha2->get_positions().get_N_sites());
    // there's no way to assert it, but we also assume they have precisely the
    // same orbitals too.  In fact, it might be useful to make a function that
    // asserts two wave functions are identical except for the particle
    // positions ...

    BOOST_ASSERT(swapped_system->subsystem_particle_counts_match());
#endif

    swapped_system->initialize(*phialpha1, *phialpha2);
}

probability_t BaseSwapPossibleWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(!transition_in_progress);
    transition_in_progress = true;

    const PositionArguments &r1 = phialpha1->get_positions();
    const PositionArguments &r2 = phialpha2->get_positions();

    chosen_particle1 = choose_random_particle(r1, rng);
    chosen_particle2 = choose_random_particle(r2, rng);

    const unsigned int particle1_destination = plan_particle_move_to_nearby_empty_site(chosen_particle1, r1, phialpha1->get_lattice(), rng);
    const unsigned int particle2_destination = plan_particle_move_to_nearby_empty_site(chosen_particle2, r2, phialpha2->get_lattice(), rng);

    // automatic reject if the subsystems will now have different particle counts of any species
    const unsigned int delta1 = calculate_subsystem_particle_change(swapped_system->get_subsystem(), r1[chosen_particle1], particle1_destination, phialpha1->get_lattice());
    const unsigned int delta2 = calculate_subsystem_particle_change(swapped_system->get_subsystem(), r2[chosen_particle2], particle2_destination, phialpha2->get_lattice());
    if (chosen_particle1.species == chosen_particle2.species) {
        if (delta1 != delta2)
            return 0;
    } else {
        if (delta1 != 0 || delta2 != 0)
            return 0;
    }

    // remember old phialpha's
    const amplitude_t old_phialpha1_psi = phialpha1->psi();
    const amplitude_t old_phialpha2_psi = phialpha2->psi();
    // implement copy-on-write
    if (!phialpha1.unique())
        phialpha1 = phialpha1->clone();
    if (!phialpha2.unique())
        phialpha2 = phialpha2->clone();
    // update phialpha's
    phialpha1->move_particle(chosen_particle1, particle1_destination);
    phialpha2->move_particle(chosen_particle2, particle2_destination);
    // determine probability ratios
    const amplitude_t phialpha1_ratio = phialpha1->psi() / old_phialpha1_psi;
    const amplitude_t phialpha2_ratio = phialpha2->psi() / old_phialpha2_psi;

    amplitude_t phibeta1_ratio(0), phibeta2_ratio(0);

    if (update_swapped_system_before_accepting) {
        // remember old phibeta's
        const amplitude_t old_phibeta1_psi = swapped_system->get_phibeta1().psi();
        const amplitude_t old_phibeta2_psi = swapped_system->get_phibeta2().psi();
        // implement copy-on-write
        if (!swapped_system.unique())
            swapped_system = boost::make_shared<SwappedSystem>(*swapped_system);
        // update phibeta's
        swapped_system->update(&chosen_particle1, &chosen_particle2, *phialpha1, *phialpha2);
        // determine probability ratios
        phibeta1_ratio = swapped_system->get_phibeta1().psi() / old_phibeta1_psi;
        phibeta2_ratio = swapped_system->get_phibeta2().psi() / old_phibeta2_psi;
    }

    // return a probability
    return probability_ratio(phialpha1_ratio, phialpha2_ratio, phibeta1_ratio, phibeta2_ratio);
}

void BaseSwapPossibleWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);
    transition_in_progress = false;

    if (!update_swapped_system_before_accepting) {
        // implement copy-on-write
        if (!swapped_system.unique())
            swapped_system = boost::make_shared<SwappedSystem>(*swapped_system);
        // update phibeta's
        swapped_system->update(&chosen_particle1, &chosen_particle2, *phialpha1, *phialpha2);
    }

#if defined(DEBUG_VMC_BASE_SWAP_POSSIBLE_WALK) || defined(DEBUG_VMC_ALL)
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
    std::cerr << std::endl;
#endif

    BOOST_ASSERT(phialpha1.unique());
    BOOST_ASSERT(phialpha2.unique());
    phialpha1->finish_particle_moved_update();
    phialpha2->finish_particle_moved_update();

    BOOST_ASSERT(swapped_system.unique());
    swapped_system->finish_update(*phialpha1, *phialpha2);
}
