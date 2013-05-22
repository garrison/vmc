#ifndef _VMC_RENYI_MOD_POSSIBLE_WALK_HPP
#define _VMC_RENYI_MOD_POSSIBLE_WALK_HPP

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
    RenyiModPossibleWalk (const boost::shared_ptr<Wavefunction<amplitude_t>::Amplitude> &wf, const boost::shared_ptr<Wavefunction<amplitude_t>::Amplitude> &wf_copy, const boost::shared_ptr<const Subsystem> &subsystem)
        : BaseSwapPossibleWalk(wf, wf_copy, subsystem, false)
        {
        }

private:
    virtual probability_t probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const override;
};

#endif
