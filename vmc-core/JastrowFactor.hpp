#ifndef _VMC_JASTROW_FACTOR_HPP
#define _VMC_JASTROW_FACTOR_HPP

#include "PositionArguments.hpp"
#include "Big.hpp"
#include "vmc-typedefs.hpp"

/**
 * Abstract base class for a projector or Jastrow factor
 */
class JastrowFactor
{
public:
    virtual ~JastrowFactor (void)
        {
        }

    // this is always going to be real and positive, but we use amplitude_t so
    // we can multiply it by other Big<amplitude_t>'s.  There is probably a
    // better way to do this.
    virtual Big<amplitude_t> compute_jastrow (const PositionArguments &r) const = 0;
};

#endif
