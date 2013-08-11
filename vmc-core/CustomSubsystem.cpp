#include <cassert>

#include "CustomSubsystem.hpp"
#include "Lattice.hpp"

bool CustomSubsystem::position_is_within (unsigned int site_index, const Lattice &lattice) const
{
    assert(lattice_makes_sense(lattice));
    assert(site_index < lattice.size());

    return this->site_status[site_index];
}

bool CustomSubsystem::lattice_makes_sense (const Lattice &lattice) const
{
    return (lattice.size() == this->site_status.size());
}
