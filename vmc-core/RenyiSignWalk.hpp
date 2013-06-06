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
    RenyiSignWalk (const std::shared_ptr<Wavefunction<amplitude_t>::Amplitude> &wf, const std::shared_ptr<Wavefunction<amplitude_t>::Amplitude> &wf_copy, const std::shared_ptr<const Subsystem> &subsystem)
        : BaseSwapPossibleWalk(wf, wf_copy, subsystem)
        {
        }

private:
    virtual probability_t probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const override;
};

#endif
