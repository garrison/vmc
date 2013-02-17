#ifndef _VMC_RENYI_SIGN_WALK_HPP
#define _VMC_RENYI_SIGN_WALK_HPP

#include "BaseSwapPossibleWalk.hpp"

/**
 * Renyi "sign" walk
 *
 * See Y. Zhang et. al., PRL 107, 067202 (2011) for explanation
 *
 * @see RenyiSignMeasurement
 */
class RenyiSignWalk : public BaseSwapPossibleWalk
{
public:
    RenyiSignWalk (const boost::shared_ptr<Wavefunction::Amplitude> &wf, const boost::shared_ptr<Wavefunction::Amplitude> &wf_copy, const boost::shared_ptr<const Subsystem> &subsystem)
        : BaseSwapPossibleWalk(wf, wf_copy, subsystem)
        {
        }

private:
    probability_t probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const;
};

#endif
