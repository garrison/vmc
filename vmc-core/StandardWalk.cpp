#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-move.hpp"

StandardWalk::StandardWalk (boost::shared_ptr<WavefunctionAmplitude> &wf_)
    : wf(wf_)
#ifndef BOOST_DISABLE_ASSERTS
    , transition_in_progress(false)
#endif
{
}

probability_t StandardWalk::compute_probability_ratio_of_random_transition (rng_class &rng)
{
    BOOST_ASSERT(!transition_in_progress);

    // remember old amplitude so we can later compute a new:old ratio
    amplitude_t old_amplitude = wf->psi();

    // move a particle and update things
    const PositionArguments &r = wf->get_positions();
    unsigned int chosen_particle = choose_random_particle(r, rng, 0);
    unsigned int new_site_index = plan_particle_move_to_nearby_empty_site(chosen_particle, r, wf->get_lattice(), rng);
    if (new_site_index == r[chosen_particle])
        return 0; // we aren't moving anything, so just reject this move
    if (!wf.unique()) // implement copy on write
        wf = wf->clone();
    wf->move_particle(chosen_particle, new_site_index);

    // calculate and return a probability
    probability_t rv = std::norm(wf->psi() / old_amplitude);
#if defined(DEBUG_VMC_STANDARD_WALK) || defined(DEBUG_VMC_ALL)
    std::cerr << "ratio " << rv << std::endl;
#endif
#ifndef BOOST_DISABLE_ASSERTS
    transition_in_progress = true;
#endif
    return rv;
}

void StandardWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    BOOST_ASSERT(wf.unique()); // ensure copy-on-write is implemented correctly
    wf->finish_particle_moved_update();

#ifndef BOOST_DISABLE_ASSERTS
    transition_in_progress = false;
#endif
}

