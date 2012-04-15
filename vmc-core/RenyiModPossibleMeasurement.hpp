#ifndef _RENYI_MOD_POSSIBLE_MEASUREMENT_HPP
#define _RENYI_MOD_POSSIBLE_MEASUREMENT_HPP

#include <cmath>

#include "Measurement.hpp"
#include "RenyiModPossibleWalk.hpp"

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
    typedef real_t measurement_value_t;

    RenyiModPossibleMeasurement (void)
        : Measurement<RenyiModPossibleWalk>(1), // we must measure after every step
          accum(0)
        {
        }

    /**
     * Returns the current value of the measurement
     */
    measurement_value_t get (void) const
        {
            return accum / (measurement_value_t) get_measurements_completed();
        }

private:
    void measure_ (const RenyiModPossibleWalk &walk)
        {
            accum += std::abs(walk.get_phibeta1().psi()
                              / walk.get_phialpha1().psi()
                              * walk.get_phibeta2().psi()
                              / walk.get_phialpha2().psi());
        }

    accumulator_t accum;
};

#endif
