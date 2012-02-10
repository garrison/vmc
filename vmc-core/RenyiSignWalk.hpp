#ifndef _RENYI_SIGN_WALK_HPP
#define _RENYI_SIGN_WALK_HPP

#include <vector>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "vmc-typedefs.hpp"
#include "Subsystem.hpp"
#include "SwappedSystem.hpp"
#include "WavefunctionAmplitude.hpp"

/**
 * Renyi "sign" walk
 *
 * See Y. Zhang et. al., PRL 107, 067202 (2011) for explanation
 *
 * @see RenyiSignMeasurement
 */
class RenyiSignWalk
{
public:
    /**
     * Constructor.
     *
     * It is essential that both wf and wf_copy have the same number of
     * particles in the subsystem.
     */
    RenyiSignWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const boost::shared_ptr<WavefunctionAmplitude> &wf_copy, boost::shared_ptr<const Subsystem> subsystem);
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

    const WavefunctionAmplitude & get_phibeta1 (void) const
        {
            return swapped_system->get_phibeta1();
        }

    const WavefunctionAmplitude & get_phibeta2 (void) const
        {
            return swapped_system->get_phibeta2();
        }

    unsigned int get_N_subsystem1 (void) const
        {
            return swapped_system->get_N_subsystem1();
        }

    unsigned int get_N_subsystem2 (void) const
        {
            return swapped_system->get_N_subsystem2();
        }

private:
    boost::shared_ptr<WavefunctionAmplitude> phialpha1, phialpha2;
    boost::shared_ptr<SwappedSystem> swapped_system;
    bool transition_in_progress;
};

#endif
