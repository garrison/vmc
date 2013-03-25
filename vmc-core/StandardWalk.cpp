#include "StandardWalk.hpp"
#include "Wavefunction.hpp"
#include "PositionArguments.hpp"

StandardWalk::StandardWalk (boost::shared_ptr<Wavefunction::Amplitude> &wfa_)
    : wfa(wfa_),
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
    const Big<amplitude_t> old_amplitude(wfa->psi());

    // choose a move and update things
    const Move move(wfa->propose_move(rng));
    if (move.size() == 0) {
        // we aren't actually moving anything, so just auto-reject this move to
        // preserve balance
        autoreject_in_progress = true;
        return 0;
    }
    if (!wfa.unique()) // implement copy on write
        wfa = wfa->clone();
    wfa->perform_move(move);

    // calculate and return a probability
    const probability_t rv = std::norm(wfa->psi().ratio(old_amplitude));
#if defined(DEBUG_VMC_STANDARD_WALK) || defined(DEBUG_VMC_ALL)
    std::cerr << "ratio " << rv << std::endl;
#endif
    return rv;
}

void StandardWalk::accept_transition (void)
{
    BOOST_ASSERT(transition_in_progress);
    BOOST_ASSERT(!autoreject_in_progress);

    BOOST_ASSERT(wfa.unique()); // ensure copy-on-write is implemented correctly
    wfa->finish_move();

    // finish_move() may recalculate the inverse from scratch, so for sanity
    // we check that the amplitude here is still nonzero.
    BOOST_ASSERT(wfa->is_nonzero());

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    transition_in_progress = false;
#endif
}

void StandardWalk::reject_transition (void)
{
    BOOST_ASSERT(transition_in_progress);

    if (!autoreject_in_progress) {
        BOOST_ASSERT(wfa.unique()); // ensure copy-on-write is implemented correctly
        wfa->cancel_move();
    }
    autoreject_in_progress = false;

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    transition_in_progress = false;
#endif
}

void StandardWalk::check_for_numerical_error (void) const
{
    wfa->check_for_numerical_error();
}
