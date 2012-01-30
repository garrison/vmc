#ifndef _RENYI_MOD_MEASUREMENT_HPP
#define _RENYI_MOD_MEASUREMENT_HPP

#include <cmath>
#include <memory>

#include "Measurement.hpp"
#include "RenyiModWalk.hpp"
#include "SwappedSystem.hpp"

/**
 * Renyi "mod" measurement
 *
 * See Y. Zhang et. al., PRL 107, 067202 (2011) for explanation.
 *
 * This measurement uses Ceperley's update method if steps_per_measurement is
 * equal to one.  Otherwise, it recalculates inverses from scratch each time.
 *
 * @see RenyiModWalk
 */
class RenyiModMeasurement : public Measurement<RenyiModWalk>
{
public:
    typedef double measurement_value_t;

    RenyiModMeasurement (const boost::shared_ptr<const Subsystem> &subsystem_, unsigned int steps_per_measurement)
        : Measurement<RenyiModWalk>(steps_per_measurement),
          subsystem(subsystem_),
          accum(0)
        {
        }

    /**
     * Returns the current value of the measurement
     */
    measurement_value_t get (void) const
        {
            return static_cast<measurement_value_t>(accum) / get_measurements_completed();
        }

private:
    void initialize_ (const RenyiModWalk &walk)
        {
            if (get_steps_per_measurement() == 1) {
                swapped_system.reset(new SwappedSystem(subsystem));
                swapped_system->initialize(walk.get_phialpha1(), walk.get_phialpha2());
            }
        }

    void measure_ (const RenyiModWalk &walk)
        {
            if (get_steps_per_measurement() == 1) {
                // update swapped_system
                std::pair<int, int> args = walk.get_swapped_system_update_args();
                swapped_system->update(args.first, args.second, walk.get_phialpha1(), walk.get_phialpha2());
                swapped_system->finish_update(walk.get_phialpha1(), walk.get_phialpha2());
                perform_the_measurement(walk);
            } else {
                swapped_system.reset(new SwappedSystem(subsystem));
                swapped_system->initialize(walk.get_phialpha1(), walk.get_phialpha2());
                perform_the_measurement(walk);
                swapped_system.reset();
            }
        }

    void repeat_measurement_ (const RenyiModWalk &walk)
        {
            if (get_steps_per_measurement() == 1) {
                perform_the_measurement(walk);
            } else {
                measure_(walk);
            }
        }

    void perform_the_measurement (const RenyiModWalk &walk)
        {
            if (swapped_system->get_N_subsystem1() == swapped_system->get_N_subsystem2()) {
                accum += std::abs(swapped_system->get_phibeta1().psi()
                                  / walk.get_phialpha1().psi()
                                  * swapped_system->get_phibeta2().psi()
                                  / walk.get_phialpha2().psi());
            }
        }

    const boost::shared_ptr<const Subsystem> subsystem;
    std::auto_ptr<SwappedSystem> swapped_system;
    accumulator_t accum;
};

#endif
