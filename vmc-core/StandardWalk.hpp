#ifndef _STANDARD_WALK_HPP
#define _STANDARD_WALK_HPP

#include <boost/shared_ptr.hpp>

#include "vmc-typedefs.hpp"
#include "Wavefunction.hpp"
#include "Walk.hpp"

/**
 * The standard walk used in VMC calculations
 *
 * (proportional to the modulus squared of the wavefunction)
 */
class StandardWalk : public Walk
{
public:
    /**
     * Constructor
     *
     * @param wf_ initial wavefunction
     */
    StandardWalk (boost::shared_ptr<Wavefunction::Amplitude> &wf_);

    /**
     * Returns the current wavefunction
     */
    const Wavefunction::Amplitude & get_wavefunction (void) const
        {
            return *wf;
        }

private:
    /**
     * Attempt a transition and return its probability ratio
     */
    probability_t compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng);

    /**
     * Accept the transition, and get the walk object into a state such that
     * new transitions can be attempted.
     */
    void accept_transition (void);

    /**
     * Reject the transition, and get the walk object into a state such that
     * new transitions can be attempted.
     */
    void reject_transition (void);

    boost::shared_ptr<Wavefunction::Amplitude> wf; // treat this as copy on write
    bool autoreject_in_progress;

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    bool transition_in_progress;
#endif
};

#endif
