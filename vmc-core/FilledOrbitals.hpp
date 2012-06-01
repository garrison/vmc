#ifndef _FILLED_ORBITALS_HPP
#define _FILLED_ORBITALS_HPP

#include <vector>

#include "OrbitalDefinitions.hpp"
#include "allowed-momentum.hpp"

/**
 * OrbitalDefinitions based on filled momenta and some given boundary conditions
 */
class FilledOrbitals : public OrbitalDefinitions
{
public:
    /**
     * Constructor
     *
     * @param momentum_sites represents which orbitals are filled
     *
     * @param lattice_ the lattice
     *
     * @param bcs the boundary conditions
     */
    FilledOrbitals (const std::vector<lw_vector<int, MAX_DIMENSION> > &momentum_sites, const boost::shared_ptr<const Lattice> &lattice_, const BoundaryConditions &bcs);

    /**
     * The boundary conditions with which the object was initialized
     */
    const BoundaryConditions boundary_conditions;
};

#endif
