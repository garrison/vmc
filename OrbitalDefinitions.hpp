#ifndef _ORBITAL_DEFINITIONS_HPP
#define _ORBITAL_DEFINITIONS_HPP

#include <boost/assert.hpp>

#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

class OrbitalDefinitions
// abstract base class
{
public:
    virtual ~OrbitalDefinitions (void)
        {
        }

    amplitude_t phi (unsigned int orbital, unsigned int site_index, const Lattice &lattice) const
        {
            BOOST_ASSERT(N_filled_orbitals < lattice.total_sites());
            BOOST_ASSERT(orbital < N_filled_orbitals);
            BOOST_ASSERT(site_index < lattice.total_sites());
            return calculate_phi(orbital, site_index, lattice);
        }

    virtual bool lattice_makes_sense (const Lattice &lattice) const = 0;

private:
    virtual amplitude_t calculate_phi (unsigned int orbital, unsigned int site_index, const Lattice &lattice) const = 0;

protected:
    OrbitalDefinitions (unsigned int N_filled_orbitals_)
        : N_filled_orbitals(N_filled_orbitals_)
        {
        }

    const unsigned int N_filled_orbitals;
};

#endif
