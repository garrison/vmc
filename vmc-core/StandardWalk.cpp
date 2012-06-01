#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-move.hpp"

StandardWalk::StandardWalk (boost::shared_ptr<WavefunctionAmplitude> &wf_)
    : wf(wf_),
      autoreject_in_progress(false)
#ifndef BOOST_DISABLE_ASSERTS
    , transition_in_progress(false)
#endif
{
}

probability_t StandardWalk::compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng)
{
    BOOST_ASSERT(!transition_in_progress);

#ifndef BOOST_DISABLE_ASSERTS
    transition_in_progress = true;
#endif

    // remember old amplitude so we can later compute a new:old ratio
    amplitude_t old_amplitude = wf->psi();

    // move a particle and update things
    const PositionArguments &r = wf->get_positions();
    Particle chosen_particle = choose_random_particle(r, rng);
    unsigned int new_site_index = plan_particle_move_to_nearby_empty_site(chosen_particle, r, wf->get_lattice(), rng);
    if (new_site_index == r[chosen_particle]) {
        // we aren't actually moving anything, so just auto-reject this move to
        // preserve balance
        autoreject_in_progress = true;
        return 0;
    }
    if (!wf.unique()) // implement copy on write
        wf = wf->clone();
    wf->perform_move(chosen_particle, new_site_index);

    // calculate and return a probability
    probability_t rv = std::norm(wf->psi() / old_amplitude);
#if defined(DEBUG_VMC_STANDARD_WALK) || defined(DEBUG_VMC_ALL)
    std::cerr << "ratio " << rv << std::endl;
#endif
    return rv;
}

void StandardWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    BOOST_ASSERT(wf.unique()); // ensure copy-on-write is implemented correctly
    if (!autoreject_in_progress)
        wf->finish_move();
    autoreject_in_progress = false;

#ifndef BOOST_DISABLE_ASSERTS
    transition_in_progress = false;
#endif
}

void StandardWalk::reject_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    BOOST_ASSERT(wf.unique()); // ensure copy-on-write is implemented correctly
    if (!autoreject_in_progress)
        wf->cancel_move();
    autoreject_in_progress = false;

#ifndef BOOST_DISABLE_ASSERTS
    transition_in_progress = false;
#endif
}
