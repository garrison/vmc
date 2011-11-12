#ifndef _RENYI_SIGN_MEASUREMENT_HPP
#define _RENYI_SIGN_MEASUREMENT_HPP

#include <cmath>

#include "RenyiSignWalk.hpp"

class RenyiSignMeasurement
{
    // fixme: in the case of a real wave function, we know that the sign will
    // always evaluate to 1 or -1.  In this case, it may make more sense to
    // store as integers how many measurements were of each value, instead of
    // using a complex<double> as an accumulator.

public:
    typedef amplitude_t measurement_value_t;

    RenyiSignMeasurement (const RenyiSignWalk &walk)
	: accum(0)
	{
	    (void) walk; // silence warning about unused variable
	}

    void measure (const RenyiSignWalk &walk)
	{
	    // we take the argument of each determinant separately instead of
	    // multiplying the determinants together first.  this is necessary
	    // because the determinants tend to be quite large, and multiplying
	    // them can lead to overflow
	    double a = 0;
	    a -= std::arg(walk.get_phialpha1().psi());
	    a -= std::arg(walk.get_phialpha2().psi());
	    a += std::arg(walk.get_phibeta1().psi());
	    a += std::arg(walk.get_phibeta2().psi());
	    const std::complex<double> i(0, 1);
#if 0
	    std::cerr << std::real(std::exp(i * a)) << "   " << std::arg(walk.get_phialpha1().psi())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phialpha2().psi())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phibeta1().psi())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phibeta2().psi())/boost::math::constants::pi<double>() << ' ' << std::endl;
#endif
	    accum += std::exp(i * a);
	}

    measurement_value_t get (unsigned int measurements_completed) const
	{
	    return accum / (measurement_value_t) measurements_completed;
	}

private:
    measurement_value_t accum; // fixme: complex accumulator_t needed

    // disable default constructor
    RenyiSignMeasurement (void);
};

#endif
