#ifndef _VMC_WALK_HPP
#define _VMC_WALK_HPP

#include <boost/noncopyable.hpp>

#include "vmc-typedefs.hpp"

class RandomNumberGenerator;

class Walk : boost::noncopyable
/** Abstract base class for a Monte Carlo walk
 *
 * @see StandardWalk for an example walk
 */
{
public:
    virtual ~Walk (void)
        {
        }

    /**
     * Attempt a transition and return its probability ratio
     */
    virtual probability_t compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng) = 0;

    /**
     * Accept the transition, and get the walk object into a state such that
     * new transitions can be attempted.
     */
    virtual void accept_transition (void) = 0;

    /**
     * Reject the transition, and get the walk object into a state such that
     * new transitions can be attempted.
     */
    virtual void reject_transition (void) = 0;

    /**
     * Performs any possible checks to see whether the Walk currently has
     * significant numerical error, and raises an exception if it does.
     *
     * This will never be called while a transition is in progress.
     */
    virtual void check_for_numerical_error (void) const = 0;
};

#endif
