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
    RenyiModPossibleWalk (std::unique_ptr<Wavefunction<amplitude_t>::Amplitude> wfa1, std::unique_ptr<Wavefunction<amplitude_t>::Amplitude> wfa2, const std::shared_ptr<const Subsystem> &subsystem)
        : BaseSwapPossibleWalk(std::move(wfa1), std::move(wfa2), subsystem, false)
        {
        }

    // CYTHON-LIMITATION: this alternative constructor exists due to poor unique_ptr support
    RenyiModPossibleWalk (std::unique_ptr<Wavefunction<amplitude_t>::Amplitude> wfa, const std::shared_ptr<const Subsystem> &subsystem)
        : BaseSwapPossibleWalk(wfa->clone(), wfa->clone(), subsystem, false)
        {
        }

private:
    virtual probability_t probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const override;
};

#endif
