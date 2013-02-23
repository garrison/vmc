#ifndef _VMC_RENYI_MOD_POSSIBLE_MEASUREMENT_HPP
#define _VMC_RENYI_MOD_POSSIBLE_MEASUREMENT_HPP

#include <cmath>

#include "Measurement.hpp"
#include "RenyiModPossibleWalk.hpp"
#include "BlockedEstimate.hpp"

/**
 * Renyi "mod/possible" measurement
 *
 * fixme: explanation needed
 *
 * @see RenyiModPossibleWalk
 */
class RenyiModPossibleMeasurement : public Measurement<RenyiModPossibleWalk>
{
public:
    RenyiModPossibleMeasurement (void)
        : Measurement<RenyiModPossibleWalk>(1) // we must measure after every step
        {
        }

    /**
     * Returns the current estimate of the measurement
     */
    const BlockedEstimate<real_t> & get_estimate (void) const
        {
            return estimate;
        }

private:
    void measure_ (const RenyiModPossibleWalk &walk)
        {
            estimate.add_value(std::abs((walk.get_phibeta1().psi()
                                         / walk.get_phialpha1().psi())
                                        * (walk.get_phibeta2().psi()
                                           / walk.get_phialpha2().psi())));
        }

    void reset (void)
        {
            estimate.reset();
        }

    BlockedEstimate<real_t> estimate;
};

#endif
