#ifndef _VMC_JASTROW_FACTOR_HPP
#define _VMC_JASTROW_FACTOR_HPP

#include "PositionArguments.hpp"
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

    virtual real_t compute_jastrow (const PositionArguments &r) const = 0;
};

#endif
