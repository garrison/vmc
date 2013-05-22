#ifndef _VMC_JASTROW_FACTOR_HPP
#define _VMC_JASTROW_FACTOR_HPP

#include "PositionArguments.hpp"
#include "Big.hpp"
#include "vmc-typedefs.hpp"

/**
 * Abstract base class for a projector or Jastrow factor
 */
template <typename AmplitudeType>
class JastrowFactor
{
public:
    virtual ~JastrowFactor (void)
        {
        }

    // this is always going to be real and positive, but we use AmplitudeType so
    // we can multiply it by other Big<AmplitudeType>'s.  There is probably a
    // better way to do this.
    virtual Big<AmplitudeType> compute_jastrow (const PositionArguments &r) const = 0;
};

#endif
