#ifndef _RENYI_MOD_POSSIBLE_WALK_HPP
#define _RENYI_MOD_POSSIBLE_WALK_HPP

#include "BaseSwapPossibleWalk.hpp"

/**
 * Renyi "mod/possible" walk
 *
 * fixme: explanation needed
 *
 * @see RenyiModPossibleMeasurement
 */
class RenyiModPossibleWalk : public BaseSwapPossibleWalk
{
public:
    RenyiModPossibleWalk (const boost::shared_ptr<Wavefunction::Amplitude> &wf, const boost::shared_ptr<Wavefunction::Amplitude> &wf_copy, boost::shared_ptr<const Subsystem> subsystem)
        : BaseSwapPossibleWalk(wf, wf_copy, subsystem, false)
        {
        }

private:
    probability_t probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const;
};

#endif
