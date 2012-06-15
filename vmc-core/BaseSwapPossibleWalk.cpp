#include <cmath>
#include <vector>

#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

#include "RandomNumberGenerator.hpp"
#include "BaseSwapPossibleWalk.hpp"
#include "random-move.hpp"

BaseSwapPossibleWalk::BaseSwapPossibleWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const boost::shared_ptr<WavefunctionAmplitude> &wf_copy, boost::shared_ptr<const Subsystem> subsystem, bool update_swapped_system_before_accepting_)
    : phialpha1(wf),
      phialpha2(wf_copy),
      swapped_system(new SwappedSystem(subsystem)),
      N_subsystem_sites(count_subsystem_sites(*subsystem, phialpha1->get_lattice())),
      update_swapped_system_before_accepting(update_swapped_system_before_accepting_),
      autoreject_in_progress(false),
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
#endif

    swapped_system->initialize(*phialpha1, *phialpha2);
}

unsigned int BaseSwapPossibleWalk::count_subsystem_sites (const Subsystem &subsystem, const Lattice &lattice)
{
    // (this is a static method)
    unsigned int N_subsystem_sites = 0;
    for (unsigned int i = 0; i < lattice.total_sites(); ++i) {
        if (subsystem.position_is_within(i, lattice))
            ++N_subsystem_sites;
    }
    return N_subsystem_sites;
}

probability_t BaseSwapPossibleWalk::compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng)
{
    BOOST_ASSERT(!transition_in_progress);
    transition_in_progress = true;

    const Lattice &lattice = phialpha1->get_lattice();
    const Subsystem &subsystem = swapped_system->get_subsystem();

    // decide which (zero-indexed) copy of the system to base this move around.
    // We call the chosen copy "copy A."
    const unsigned int copy_A = rng.random_small_uint(2);

    const PositionArguments &r_A = (copy_A == 0) ? phialpha1->get_positions() : phialpha2->get_positions();
    const PositionArguments &r_B = (copy_A == 0) ? phialpha2->get_positions() : phialpha1->get_positions();

    // first choose a random move in copy A
    chosen_particle_A = choose_random_particle(r_A, rng);
    const unsigned int particle_A_destination = plan_particle_move_to_nearby_empty_site(chosen_particle_A, r_A, lattice, rng);

    if (r_A[chosen_particle_A] == particle_A_destination) {
        autoreject_in_progress = true;
        return 0;
    }

    // If the particle we are moving in copy A is not going to change its
    // subsystem status, then we go ahead and make that move.  If, however, the
    // particle we are moving in copy A changes its subsystem status (i.e.,
    // enters or leaves the subsystem), then we must also choose a move in copy
    // B that does likewise, and take into account the transition probability
    // ratio so that the simulation maintains balance.

    unsigned int particle_B_destination = -1; // uninitialized
    real_t transition_ratio(1);

    const int copy_A_subsystem_particle_change = calculate_subsystem_particle_change(swapped_system->get_subsystem(), r_A[chosen_particle_A], particle_A_destination, lattice);
    const unsigned int species = chosen_particle_A.species;
    if (copy_A_subsystem_particle_change != 0) {
        const bool candidate_particle_B_subsystem_status = (copy_A_subsystem_particle_change == -1);

        // determine reverse/forward transition attempt probability ratios
        const unsigned int N_sites = r_A.get_N_sites();
        const unsigned int N_filled = r_A.get_N_filled(species);
        const unsigned int N_within_subsystem = swapped_system->get_N_subsystem(species);
        const unsigned int N_outside_subsystem = N_filled - N_within_subsystem;
        const unsigned int N_vacant_within_subsystem = N_subsystem_sites - N_within_subsystem;
        const unsigned int N_vacant_outside_subsystem = (N_sites - N_filled) - N_vacant_within_subsystem;
        unsigned int forward_particle_possibilities, forward_vacant_possibilities,
            reverse_particle_possibilities, reverse_vacant_possibilities;
        if (copy_A_subsystem_particle_change == 1) {
            forward_particle_possibilities = N_outside_subsystem;
            forward_vacant_possibilities = N_vacant_within_subsystem;
            reverse_particle_possibilities = N_within_subsystem + 1;
            reverse_vacant_possibilities = N_vacant_outside_subsystem + 1;
        } else {
            BOOST_ASSERT(copy_A_subsystem_particle_change == -1);
            forward_particle_possibilities = N_within_subsystem;
            forward_vacant_possibilities = N_vacant_outside_subsystem;
            reverse_particle_possibilities = N_outside_subsystem + 1;
            reverse_vacant_possibilities = N_vacant_within_subsystem + 1;
        }
        if (forward_particle_possibilities == 0 || forward_vacant_possibilities == 0) {
            autoreject_in_progress = true;
            return 0;
        }
        transition_ratio = real_t(forward_particle_possibilities * forward_vacant_possibilities) / real_t(reverse_particle_possibilities * reverse_vacant_possibilities);

        // choose a particle from B with the same subsystem status as the
        // particle we are moving in A
        std::vector<unsigned int> candidate_particle_B_array;
        candidate_particle_B_array.reserve(forward_particle_possibilities);
        for (unsigned int i = 0; i < N_filled; ++i) {
            if (subsystem.position_is_within(r_B[Particle(i, species)], lattice) == candidate_particle_B_subsystem_status)
                candidate_particle_B_array.push_back(i);
        }
        BOOST_ASSERT(forward_particle_possibilities == candidate_particle_B_array.size());
        chosen_particle_B = Particle(candidate_particle_B_array[rng.random_small_uint(candidate_particle_B_array.size())], species);

        // choose a destination such that the particle in copy B will change its
        // subsystem status
        const bool destination_B_in_subsystem = !candidate_particle_B_subsystem_status;
        std::vector<unsigned int> candidate_destination_B_array;
        candidate_destination_B_array.reserve(forward_vacant_possibilities);
        for (unsigned int i = 0; i < N_sites; ++i) {
            if (!r_B.is_occupied(i, species) && subsystem.position_is_within(i, lattice) == destination_B_in_subsystem)
            candidate_destination_B_array.push_back(i);
        }
        BOOST_ASSERT(forward_vacant_possibilities == candidate_destination_B_array.size());
        particle_B_destination = candidate_destination_B_array[rng.random_small_uint(candidate_destination_B_array.size())];
    }

    // translate from "copy A or B" language back to "copy 1 or 2"
    const Particle * const chosen_particle_B_ptr = (copy_A_subsystem_particle_change != 0) ? &chosen_particle_B : 0;
    chosen_particle1 = (copy_A == 0) ? &chosen_particle_A : chosen_particle_B_ptr;
    chosen_particle2 = (copy_A == 0) ? chosen_particle_B_ptr : &chosen_particle_A;

    // move particles, determining phialpha probability ratios
    amplitude_t phialpha1_ratio(1), phialpha2_ratio(1);
    if (chosen_particle1) {
        const amplitude_t old_phialpha1_psi = phialpha1->psi();
        if (!phialpha1.unique()) // copy-on-write
            phialpha1 = phialpha1->clone();
        Move move;
        move.push_back(SingleParticleMove(*chosen_particle1, (copy_A == 0) ? particle_A_destination : particle_B_destination));
        phialpha1->perform_move(move);
        phialpha1_ratio = phialpha1->psi() / old_phialpha1_psi;
    }
    if (chosen_particle2) {
        const amplitude_t old_phialpha2_psi = phialpha2->psi();
        if (!phialpha2.unique()) // copy-on-write
            phialpha2 = phialpha2->clone();
        Move move;
        move.push_back(SingleParticleMove(*chosen_particle2, (copy_A == 0) ? particle_B_destination : particle_A_destination));
        phialpha2->perform_move(move);
        phialpha2_ratio = phialpha2->psi() / old_phialpha2_psi;
    }

    amplitude_t phibeta1_ratio(0), phibeta2_ratio(0);

    if (update_swapped_system_before_accepting) {
        // remember old phibeta's
        const amplitude_t old_phibeta1_psi = swapped_system->get_phibeta1().psi();
        const amplitude_t old_phibeta2_psi = swapped_system->get_phibeta2().psi();
        // implement copy-on-write
        if (!swapped_system.unique())
            swapped_system = boost::make_shared<SwappedSystem>(*swapped_system);
        // update phibeta's
        swapped_system->update(chosen_particle1, chosen_particle2, *phialpha1, *phialpha2);
        // determine probability ratios
        phibeta1_ratio = swapped_system->get_phibeta1().psi() / old_phibeta1_psi;
        phibeta2_ratio = swapped_system->get_phibeta2().psi() / old_phibeta2_psi;
    }

    // return a probability
    return transition_ratio * probability_ratio(phialpha1_ratio, phialpha2_ratio, phibeta1_ratio, phibeta2_ratio);
}

void BaseSwapPossibleWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);
    transition_in_progress = false;

    BOOST_ASSERT(!autoreject_in_progress);

    if (!update_swapped_system_before_accepting) {
        // implement copy-on-write
        if (!swapped_system.unique())
            swapped_system = boost::make_shared<SwappedSystem>(*swapped_system);
        // update phibeta's
        swapped_system->update(chosen_particle1, chosen_particle2, *phialpha1, *phialpha2);
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

    if (chosen_particle1) {
        BOOST_ASSERT(phialpha1.unique());
        phialpha1->finish_move();
    }
    if (chosen_particle2) {
        BOOST_ASSERT(phialpha2.unique());
        phialpha2->finish_move();
    }

    BOOST_ASSERT(swapped_system.unique());
    swapped_system->finish_update(*phialpha1, *phialpha2);
}

void BaseSwapPossibleWalk::reject_transition (void)
{
    BOOST_ASSERT(transition_in_progress);
    transition_in_progress = false;

    if (autoreject_in_progress) {
        autoreject_in_progress = false;
        return;
    }

    // as usual, we modify (restore, in this case) the phialpha's before
    // telling the swapped_system to restore itself
    if (chosen_particle1) {
        BOOST_ASSERT(phialpha1.unique());
        phialpha1->cancel_move();
    }
    if (chosen_particle2) {
        BOOST_ASSERT(phialpha2.unique());
        phialpha2->cancel_move();
    }

    if (update_swapped_system_before_accepting) {
        // ensure copy-on-write was implemented
        BOOST_ASSERT(swapped_system.unique());
        // cancel the phibeta updates
        swapped_system->cancel_update(*phialpha1, *phialpha2);
    }
}
