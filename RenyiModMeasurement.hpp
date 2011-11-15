#ifndef _RENYI_MOD_MEASUREMENT_HPP
#define _RENYI_MOD_MEASUREMENT_HPP

#include <cmath>
#include <memory>

#include "Measurement.hpp"
#include "RenyiModWalk.hpp"
#include "SwappedSystem.hpp"

class RenyiModMeasurement : public Measurement<RenyiModWalk>
{
public:
    typedef double measurement_value_t;

    RenyiModMeasurement (const boost::shared_ptr<const Subsystem> &subsystem)
	: swapped_system(new SwappedSystem(subsystem)),
	  accum(0)
	{
	}

    measurement_value_t get (void) const
	{
	    return static_cast<measurement_value_t>(accum) / get_measurements_completed();
	}

private:
    void initialize_ (const RenyiModWalk &walk)
	{
	    swapped_system->initialize(walk.get_phialpha1(), walk.get_phialpha2());
	}

    void measure_ (const RenyiModWalk &walk)
	{
	    // NOTE: this must be called on *every* accepted step of RenyiModWalk.
	    // Otherwise the phibeta's will get messed up.

	    // update swapped_system
	    std::pair<int, int> args = walk.get_swapped_system_update_args();
	    swapped_system->update(args.first, args.second, walk.get_phialpha1(), walk.get_phialpha2());
	    swapped_system->finish_update(walk.get_phialpha1(), walk.get_phialpha2());

	    perform_the_measurement(walk);
	}

    void repeat_measurement_ (const RenyiModWalk &walk)
	{
	    perform_the_measurement(walk);
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

    // disable default constructor
    RenyiModMeasurement (void);

    std::auto_ptr<SwappedSystem> swapped_system;
    accumulator_t accum;
};

#endif
