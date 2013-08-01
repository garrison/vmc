#ifndef _VMC_STANDARD_WALK_HPP
#define _VMC_STANDARD_WALK_HPP

#include <memory>

#include "vmc-typedefs.hpp"
#include "vmc-real-part.hpp"
#include "Wavefunction.hpp"
#include "Walk.hpp"

/**
 * The standard walk used in VMC calculations
 *
 * (proportional to the modulus squared of the wavefunction)
 */
template <typename AmplitudeType>
class StandardWalk : public Walk<typename RealPart<AmplitudeType>::type>
{
public:
    typedef typename RealPart<AmplitudeType>::type ProbabilityType;

    /**
     * Constructor
     *
     * @param wfa_ initial wavefunctionamplitude
     */
    StandardWalk (std::unique_ptr<typename Wavefunction<AmplitudeType>::Amplitude> wfa_);

    /**
     * Returns the wavefunctionamplitude of this walk.  The returned object
     * will be updated as the walk is performed.
     */
    const typename Wavefunction<AmplitudeType>::Amplitude & get_wavefunctionamplitude (void) const
        {
            return *wfa;
        }

private:
    /**
     * Attempt a transition and return its probability ratio
     */
    virtual ProbabilityType compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng) override;

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

    std::unique_ptr<typename Wavefunction<AmplitudeType>::Amplitude> wfa;
    bool autoreject_in_progress;

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    bool transition_in_progress;
#endif
};

#include "vmc-real-part.hpp"
#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class StandardWalk<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif
