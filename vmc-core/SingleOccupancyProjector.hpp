#ifndef _SINGLE_OCCUPANCY_PROJECTOR_HPP
#define _SINGLE_OCCUPANCY_PROJECTOR_HPP

#include "JastrowFactor.hpp"

/**
 * Projects out doubly-occupied sites.
 *
 * Currently, should only be used when get_N_species == 2.
 */
class SingleOccupancyProjector : public JastrowFactor
{
    real_t compute_jastrow (const PositionArguments &r) const;
};

#endif
