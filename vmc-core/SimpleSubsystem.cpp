#include <cassert>

#include "SimpleSubsystem.hpp"

bool SimpleSubsystem::position_is_within (unsigned int site_index, const Lattice &lattice) const
{
    assert(lattice_makes_sense(lattice));

    const LatticeSite site(lattice[site_index]);
    for (unsigned int i = 0; i < lattice.n_dimensions(); ++i) {
        assert(site[i] >= 0);
        if (site[i] >= (int) subsystem_length[i])
            return false;
    }
    return true;
}

bool SimpleSubsystem::lattice_makes_sense (const Lattice &lattice) const
{
    if (lattice.n_dimensions() != subsystem_length.size())
        return false;
    for (unsigned int i = 0; i < lattice.n_dimensions(); ++i) {
        if (lattice.dimensions[i] < (int) subsystem_length[i])
            return false;
    }
    return true;
}
