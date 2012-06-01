#ifndef _SIMPLE_SUBSYSTEM_HPP
#define _SIMPLE_SUBSYSTEM_HPP

#include "lw_vector.hpp"
#include "Subsystem.hpp"
#include "Lattice.hpp"

/**
 * Represents any subsystem that is a parallelpiped aligned with the lattice's
 * primitive vectors
 */
class SimpleSubsystem : public Subsystem
{
public:
    explicit SimpleSubsystem (const lw_vector<unsigned int, MAX_DIMENSION> &subsystem_length_)
        : subsystem_length(subsystem_length_)
        {
        }

    bool position_is_within (unsigned int site_index, const Lattice &lattice) const;

    bool lattice_makes_sense (const Lattice &lattice) const;

private:
    const lw_vector<unsigned int, MAX_DIMENSION> subsystem_length;
};

#endif
