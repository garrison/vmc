#ifndef _VMC_RENYI_SIGN_MEASUREMENT_HPP
#define _VMC_RENYI_SIGN_MEASUREMENT_HPP

#include <cmath>

#include "Measurement.hpp"
#include "RenyiSignWalk.hpp"
#include "BlockedEstimate.hpp"

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
    // using amplitude_t as an accumulator.

public:
    RenyiSignMeasurement (void)
        : Measurement<RenyiSignWalk>(1) // we must measure after every step
        {
        }

    /**
     * Returns the current estimate of the measurement
     */
    const BlockedEstimate<amplitude_t> & get_estimate (void) const
        {
            return estimate;
        }

private:
    virtual void measure_ (const RenyiSignWalk &walk) override
        {
            // we take the argument of each determinant separately instead of
            // multiplying the determinants together first.  this is necessary
            // because the determinants tend to be quite large, and multiplying
            // them can lead to overflow
            using std::arg;
            using std::exp;
            real_t a = 0;
            a -= arg(walk.get_phialpha1().psi().get_base());
            a -= arg(walk.get_phialpha2().psi().get_base());
            a += arg(walk.get_phibeta1().psi().get_base());
            a += arg(walk.get_phibeta2().psi().get_base());
            const complex_t i(0, 1);
#if 0
            std::cerr << std::real(std::exp(i * a)) << "   " << std::arg(walk.get_phialpha1().psi())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phialpha2().psi())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phibeta1().psi())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phibeta2().psi())/boost::math::constants::pi<double>() << ' ' << std::endl;
#endif
            estimate.add_value(exp(i * a));
        }

    BlockedEstimate<amplitude_t> estimate;
};

#endif
