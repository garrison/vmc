#ifndef _WALK_HPP
#define _WALK_HPP

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
};

#endif
