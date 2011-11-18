#ifndef _CHAIN1D_ORBITALS_HPP
#define _CHAIN1D_ORBITALS_HPP

#include "OrbitalDefinitions.hpp"

class Chain1dOrbitals : public OrbitalDefinitions
{
public:
    Chain1dOrbitals (unsigned int N_filled)
        : OrbitalDefinitions(N_filled)
        {
        }

    bool lattice_makes_sense (const Lattice &lattice) const;

private:
    amplitude_t calculate_phi (unsigned int orbital, unsigned int site_index, const Lattice &lattice) const;
};

#endif
