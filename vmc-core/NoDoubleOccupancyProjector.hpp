#ifndef _VMC_NO_DOUBLE_OCCUPANCY_PROJECTOR_HPP
#define _VMC_NO_DOUBLE_OCCUPANCY_PROJECTOR_HPP

#include "JastrowFactor.hpp"

/**
 * Projects out doubly-occupied sites.
 *
 * Currently, should only be used when get_N_species == 2.
 */
class NoDoubleOccupancyProjector : public JastrowFactor
{
    real_t compute_jastrow (const PositionArguments &r) const;
};

#endif
