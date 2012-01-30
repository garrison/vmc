#ifndef _RENYI_SIGN_MEASUREMENT_HPP
#define _RENYI_SIGN_MEASUREMENT_HPP

#include <cmath>

#include "Measurement.hpp"
#include "RenyiSignWalk.hpp"

/**
 * Renyi "sign" measurement
 *
 * See Y. Zhang et. al., PRL 107, 067202 (2011) for explanation
 *
 * @see RenyiSignWalk
 */
class RenyiSignMeasurement : public Measurement<RenyiSignWalk>
{
    // fixme: in the case of a real wave function, we know that the sign will
    // always evaluate to 1 or -1.  In this case, it may make more sense to
    // store as integers how many measurements were of each value, instead of
    // using a complex<double> as an accumulator.

public:
    typedef amplitude_t measurement_value_t;

    RenyiSignMeasurement (void)
        : Measurement<RenyiSignWalk>(1), // we must measure after every step
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
    void measure_ (const RenyiSignWalk &walk)
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

    measurement_value_t accum; // fixme: complex accumulator_t needed
};

#endif
