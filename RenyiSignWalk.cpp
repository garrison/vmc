#include <cmath>
#include <vector>

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>
#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

#include "RenyiSignWalk.hpp"
#include "random-move.hpp"

RenyiSignWalk::RenyiSignWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const boost::shared_ptr<WavefunctionAmplitude> &wf_copy, boost::shared_ptr<const Subsystem> subsystem)
    : phialpha1(wf),
      phialpha2(wf_copy),
      swapped_system(new SwappedSystem(subsystem)),
      transition_in_progress(false)
{
    BOOST_ASSERT(&phialpha1->get_lattice() == &phialpha2->get_lattice());
    BOOST_ASSERT(phialpha1->get_positions().get_N_filled() == phialpha2->get_positions().get_N_filled());
    BOOST_ASSERT(phialpha1->get_positions().get_N_sites() == phialpha2->get_positions().get_N_sites());
    // there's no way to assert it, but we also assume they have precisely the
    // same orbitals too.  In fact, it might be useful to make a function that
    // asserts two wave functions are identical except for the particle
    // positions ...

    BOOST_ASSERT(swapped_system->get_N_subsystem1() == swapped_system->get_N_subsystem2());

    swapped_system->initialize(*phialpha1, *phialpha2);

    // count how many sites of the lattice are in the subsystem.  we will need
    // this later.
    const Lattice &lattice = phialpha1->get_lattice();
    N_subsystem_sites = 0;
    for (unsigned int i = 0; i < lattice.total_sites(); ++i) {
        if (subsystem->position_is_within(i, lattice))
            ++N_subsystem_sites;
    }
}

probability_t RenyiSignWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(!transition_in_progress);
    transition_in_progress = true;

    const Lattice &lattice = phialpha1->get_lattice();
    const Subsystem &subsystem = swapped_system->get_subsystem();

    // decide which copy of the system to base this move around.  We call the
    // chosen copy "copy A."
    boost::uniform_smallint<> copy_distribution(1, 2);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > copy_gen(rng, copy_distribution);
    unsigned int copy_A = copy_gen();

    const PositionArguments &r_A = (copy_A == 1) ? phialpha1->get_positions() : phialpha2->get_positions();
    const PositionArguments &r_B = (copy_A == 1) ? phialpha2->get_positions() : phialpha1->get_positions();

    // first choose a random move in copy A
    const unsigned int chosen_particle_A = choose_random_particle(r_A, rng);
    const unsigned int particle_A_destination = plan_particle_move_to_nearby_empty_site(chosen_particle_A, r_A, lattice, rng);

    // now we need to come up with a move in copy B that results in the same
    // change in subsystem particle number as the move in copy A.  And we need
    // to compute the transition probability ratio so the simulation maintains
    // balance.
    unsigned int chosen_particle_B, particle_B_destination;
    real_t transition_ratio;

    int copy_A_subsystem_particle_change = calculate_subsystem_particle_change(swapped_system->get_subsystem(), r_A[chosen_particle_A], particle_A_destination, lattice);
    if (copy_A_subsystem_particle_change == 0) {
        chosen_particle_B = choose_random_particle(r_B, rng);
        transition_ratio = 1;
    } else {
        const bool candidate_particle_B_subsystem_status = (copy_A_subsystem_particle_change == -1);
        std::vector<unsigned int> candidate_particle_B_array;
        candidate_particle_B_array.reserve(r_B.size());
        for (unsigned int i = 0; i < r_B.size(); ++i) {
            if (subsystem.position_is_within(r_B[i], lattice) == candidate_particle_B_subsystem_status)
                candidate_particle_B_array.push_back(i);
        }
        if (candidate_particle_B_array.empty())
            return 0;

        boost::uniform_smallint<> candidate_distribution(0, candidate_particle_B_array.size() - 1);
        boost::variate_generator<rng_class&, boost::uniform_smallint<> > candidate_gen(rng, candidate_distribution);
        chosen_particle_B = candidate_particle_B_array[candidate_gen()];

        // determine reverse/forward transition attempt probability ratio
        BOOST_ASSERT(swapped_system->get_N_subsystem1() == swapped_system->get_N_subsystem2());
        const unsigned int N_subsystem = swapped_system->get_N_subsystem1();
        const unsigned int N_filled = r_A.get_N_filled();
        const unsigned int N_sites = r_A.get_N_sites();
        const unsigned int forward_particle_possibilities = (copy_A_subsystem_particle_change == 1) ? (N_filled - N_subsystem) : N_subsystem;
        const unsigned int forward_vacant_possibilities = (copy_A_subsystem_particle_change == 1) ? (N_subsystem_sites - N_subsystem) : (N_sites - N_subsystem_sites - N_filled + N_subsystem);
        const unsigned int reverse_particle_possibilities = (copy_A_subsystem_particle_change == 1) ? (N_subsystem + 1) : (N_filled - N_subsystem + 1);
        const unsigned int reverse_vacant_possibilities = (copy_A_subsystem_particle_change == 1) ? (N_sites - N_subsystem_sites - N_filled + N_subsystem + 1) : (N_subsystem_sites - N_subsystem + 1);
        transition_ratio = real_t(forward_particle_possibilities * forward_vacant_possibilities) / real_t(reverse_particle_possibilities * reverse_vacant_possibilities);
    }

    // choose a destination such that the particle number in copy B remains
    // equal to the particle number in copy A
    bool destination_B_in_subsystem = ((copy_A_subsystem_particle_change == 1)
                                       || (copy_A_subsystem_particle_change == 0
                                           && subsystem.position_is_within(r_B[chosen_particle_B], lattice)));
    std::vector<unsigned int> candidate_destination_B_array;
    candidate_destination_B_array.reserve(r_B.get_N_sites() - r_B.get_N_filled());
    for (unsigned int i = 0; i < r_B.get_N_sites(); ++i) {
        if (!r_B.is_occupied(i) && subsystem.position_is_within(i, lattice) == destination_B_in_subsystem)
            candidate_destination_B_array.push_back(i);
    }
    if (candidate_destination_B_array.empty())
        return 0;

    boost::uniform_smallint<> destination_distribution(0, candidate_destination_B_array.size() - 1);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > destination_gen(rng, destination_distribution);
    particle_B_destination = candidate_destination_B_array[destination_gen()];

    // translate from "copy A or B" language back to "copy 1 or 2"
    const unsigned int chosen_particle1 = (copy_A == 1) ? chosen_particle_A : chosen_particle_B;
    const unsigned int chosen_particle2 = (copy_A == 1) ? chosen_particle_B : chosen_particle_A;
    const unsigned int particle1_destination = (copy_A == 1) ? particle_A_destination : particle_B_destination;
    const unsigned int particle2_destination = (copy_A == 1) ? particle_B_destination : particle_A_destination;

    // remember old wavefunction amplitudes so we can compute ratios
    const amplitude_t old_phialpha1_psi = phialpha1->psi();
    const amplitude_t old_phialpha2_psi = phialpha2->psi();
    const amplitude_t old_phibeta1_psi = swapped_system->get_phibeta1().psi();
    const amplitude_t old_phibeta2_psi = swapped_system->get_phibeta2().psi();

    // implement copy-on-write
    if (!phialpha1.unique())
        phialpha1 = phialpha1->clone();
    if (!phialpha2.unique())
        phialpha2 = phialpha2->clone();
    if (!swapped_system.unique())
        swapped_system = boost::make_shared<SwappedSystem>(*swapped_system);

    // update phialpha's
    phialpha1->move_particle(chosen_particle1, particle1_destination);
    phialpha2->move_particle(chosen_particle2, particle2_destination);

    // update phibeta's
    swapped_system->update(chosen_particle1, chosen_particle2, *phialpha1, *phialpha2);

    // return a probability
    const amplitude_t phialpha1_ratio = phialpha1->psi() / old_phialpha1_psi;
    const amplitude_t phialpha2_ratio = phialpha2->psi() / old_phialpha2_psi;
    const amplitude_t phibeta1_ratio = swapped_system->get_phibeta1().psi() / old_phibeta1_psi;
    const amplitude_t phibeta2_ratio = swapped_system->get_phibeta2().psi() / old_phibeta2_psi;
    return std::abs(phialpha1_ratio * phialpha2_ratio * phibeta1_ratio * phibeta2_ratio) * transition_ratio;
}

void RenyiSignWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);
    transition_in_progress = false;

#if defined(DEBUG_VMC_RENYI_SIGN_WALK) || defined(DEBUG_VMC_ALL)
    const PositionArguments &r1 = phialpha1->get_positions();
    for (unsigned int i = 0; i < r1.size(); ++i)
        std::cerr << r1[i] << ' ';
    std::cerr << std::endl;

    const PositionArguments &r2 = phialpha2->get_positions();
    for (unsigned int i = 0; i < r2.size(); ++i)
        std::cerr << r2[i] << ' ';
    std::cerr << std::endl << swapped_system->get_N_subsystem1() << std::endl;
#endif

    BOOST_ASSERT(phialpha1.unique());
    BOOST_ASSERT(phialpha2.unique());
    phialpha1->finish_particle_moved_update();
    phialpha2->finish_particle_moved_update();

    BOOST_ASSERT(swapped_system.unique());
    swapped_system->finish_update(*phialpha1, *phialpha2);
}
