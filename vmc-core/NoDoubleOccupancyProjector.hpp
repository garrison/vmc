#ifndef _VMC_NO_DOUBLE_OCCUPANCY_PROJECTOR_HPP
#define _VMC_NO_DOUBLE_OCCUPANCY_PROJECTOR_HPP

#include "JastrowFactor.hpp"

// FIXME: this should be implemented as a Projector, not a JastrowFactor.

/**
 * Projects out doubly-occupied sites.
 *
 * Currently, should only be used when get_N_species == 2.
 */
template <typename AmplitudeType>
class NoDoubleOccupancyProjector : public JastrowFactor<AmplitudeType>
{
    virtual Big<AmplitudeType> compute_jastrow (const PositionArguments &r) const override;
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class NoDoubleOccupancyProjector<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif
