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
    std::pair<int, int> get_swapped_system_update_args (void) const
        {
            BOOST_ASSERT(transition_copy_in_progress == 0);
            // REMEMBER: this function also assumes that accept_transition()
            // has been called at least once.
            return swapped_system_update_args;
        }

private:
    boost::shared_ptr<WavefunctionAmplitude> phialpha1, phialpha2;
    unsigned int transition_copy_in_progress;
    unsigned int chosen_particle;
    std::pair<int, int> swapped_system_update_args;
};

#endif
