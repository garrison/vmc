#include <boost/assert.hpp>

#include "NoDoubleOccupancyProjector.hpp"

real_t NoDoubleOccupancyProjector::compute_jastrow (const PositionArguments &r) const
{
    // FIXME: we should make this requirement more explicit to callers
    BOOST_ASSERT(r.get_N_species() == 2);

    for (unsigned int i = 0; i < r.get_N_filled(0); ++i) {
        if (r.is_occupied(r[Particle(i, 0)], 1))
            return 0;
    }
    return 1;
}
