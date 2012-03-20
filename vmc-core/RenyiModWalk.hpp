#ifndef _RENYI_MOD_WALK_HPP
#define _RENYI_MOD_WALK_HPP

#include <utility>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "vmc-typedefs.hpp"
#include "WavefunctionAmplitude.hpp"

/**
 * Renyi "mod" walk
 *
 * See Y. Zhang et. al., PRL 107, 067202 (2011) for explanation
 *
 * @see RenyiModMeasurement
 */
class RenyiModWalk
{
public:
    RenyiModWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const boost::shared_ptr<WavefunctionAmplitude> &wf_copy);
    probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
    void accept_transition (void);

    const WavefunctionAmplitude & get_phialpha1 (void) const
        {
            return *phialpha1;
        }

    const WavefunctionAmplitude & get_phialpha2 (void) const
        {
            return *phialpha2;
        }

    /**
     * Returns the first two arguments that should be passed to
     * SwappedSystem::update().  This should be done by RenyiModMeasurement
     * after each transition.
     */
    std::pair<const Particle *, const Particle *> get_swapped_system_update_args (void) const
        {
            // make sure transition has been accepted
            BOOST_ASSERT(transition_copy_in_progress == 0);

            BOOST_ASSERT(transition_copy_just_completed == 1
                         || transition_copy_just_completed == 2);

            if (transition_copy_just_completed == 1)
                return std::make_pair(&chosen_particle, (Particle *) 0);
            else
                return std::make_pair((Particle *) 0, &chosen_particle);
        }

private:
    boost::shared_ptr<WavefunctionAmplitude> phialpha1, phialpha2;
    int transition_copy_in_progress, transition_copy_just_completed;
    Particle chosen_particle;
};

#endif
