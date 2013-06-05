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
    virtual void measure_ (const RenyiModPossibleWalk &walk) override
        {
            using std::abs;
            estimate.add_value(abs(walk.get_phibeta1().psi().ratio(walk.get_phialpha1().psi())
                                   * walk.get_phibeta2().psi().ratio(walk.get_phialpha2().psi())));
        }

    BlockedEstimate<real_t> estimate;
};

#endif
