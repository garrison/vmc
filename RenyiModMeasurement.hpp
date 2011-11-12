#ifndef _RENYI_MOD_MEASUREMENT_HPP
#define _RENYI_MOD_MEASUREMENT_HPP

#include <cmath>

#include "RenyiModWalk.hpp"

class RenyiModMeasurement
{
public:
    typedef std::vector<double> measurement_value_t;

    RenyiModMeasurement (const RenyiModWalk &walk)
	: accum(walk.get_subsystem_array_size())
	{
	}

    void measure (const RenyiModWalk &walk)
	{
	    for (unsigned int i = 0; i < accum.size(); ++i) {
		if (walk.get_N_subsystem1(i) == walk.get_N_subsystem2(i)) {
		    accum[i] += std::abs(walk.get_phibeta1(i).psi()
					 / walk.get_phialpha1().psi()
					 * walk.get_phibeta2(i).psi()
					 / walk.get_phialpha2().psi());
		}
	    }
	}

    measurement_value_t get (unsigned int measurements_completed) const
	{
	    std::vector<double> rv(accum.size());
	    for (unsigned int i = 0; i < accum.size(); ++i)
		rv[i] = static_cast<double>(accum[i]) / measurements_completed;
	    return rv;
	}

private:
    std::vector<accumulator_t> accum;

    // disable default constructor
    RenyiModMeasurement (void);
};

#endif
