#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"

StandardWalk::StandardWalk (boost::shared_ptr<WavefunctionAmplitude> &wf_)
    : wf(wf_),
      autoreject_in_progress(false)
#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    , transition_in_progress(false)
#endif
{
}

probability_t StandardWalk::compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng)
{
    BOOST_ASSERT(!transition_in_progress);

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    transition_in_progress = true;
#endif

    // remember old amplitude so we can later compute a new:old ratio
    amplitude_t old_amplitude = wf->psi();

    // choose a move and update things
    const Move move(wf->propose_move(rng));
    if (move.size() == 0) {
        // we aren't actually moving anything, so just auto-reject this move to
        // preserve balance
        autoreject_in_progress = true;
        return 0;
    }
    if (!wf.unique()) // implement copy on write
        wf = wf->clone();
    wf->perform_move(move);

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
    BOOST_ASSERT(!autoreject_in_progress);

    BOOST_ASSERT(wf.unique()); // ensure copy-on-write is implemented correctly
    wf->finish_move();

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    transition_in_progress = false;
#endif
}

void StandardWalk::reject_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    if (!autoreject_in_progress) {
        BOOST_ASSERT(wf.unique()); // ensure copy-on-write is implemented correctly
        wf->cancel_move();
    }
    autoreject_in_progress = false;

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    transition_in_progress = false;
#endif
}
