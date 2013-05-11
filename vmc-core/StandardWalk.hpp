#ifndef _VMC_STANDARD_WALK_HPP
#define _VMC_STANDARD_WALK_HPP

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
     * @param wfa_ initial wavefunction
     */
    StandardWalk (boost::shared_ptr<Wavefunction::Amplitude> &wfa_);

    /**
     * Returns the current wavefunction
     */
    const Wavefunction::Amplitude & get_wavefunctionamplitude (void) const
        {
            return *wfa;
        }

private:
    /**
     * Attempt a transition and return its probability ratio
     */
    virtual probability_t compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng) override;

    /**
     * Accept the transition, and get the walk object into a state such that
     * new transitions can be attempted.
     */
    virtual void accept_transition (void) override;

    /**
     * Reject the transition, and get the walk object into a state such that
     * new transitions can be attempted.
     */
    virtual void reject_transition (void) override;

    virtual void check_for_numerical_error (void) const override;

    boost::shared_ptr<Wavefunction::Amplitude> wfa; // treat this as copy on write
    bool autoreject_in_progress;

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    bool transition_in_progress;
#endif
};

#endif
