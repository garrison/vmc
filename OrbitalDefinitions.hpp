#ifndef _ORBITAL_DEFINITIONS_HPP
#define _ORBITAL_DEFINITIONS_HPP

#include <Eigen/Dense>
#include <boost/assert.hpp>

#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

/**
 * "Abstract base class" for defining the orbitals in a determinantal
 * wavefunction
 *
 * (It's not strictly an abstract base class, as it has no undefined virtual
 * functions, but we should treat it as such.)
 */
class OrbitalDefinitions
{
protected:
    OrbitalDefinitions (unsigned int N_filled_orbitals, const boost::shared_ptr<const Lattice> &lattice_)
        : orbitals(N_filled_orbitals, lattice_->total_sites()),
          lattice(lattice_)
        {
        }

public:
    virtual ~OrbitalDefinitions (void)
        {
        }

    /**
     * Returns the set of orbitals at the given position
     */
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic>::ConstColXpr at_position (unsigned int site_index) const
        {
            BOOST_ASSERT(site_index < orbitals.cols());
            return orbitals.col(site_index);
        }

    unsigned int get_N_filled (void) const
        {
            return orbitals.rows();
        }

    unsigned int get_N_sites (void) const
        {
            return orbitals.cols();
        }

    const boost::shared_ptr<const Lattice> & get_lattice_ptr (void) const
        {
            return lattice;
        }

protected:
    /**
     * Stores the definition of each orbital at each possible position
     *
     * The rows of the matrix represent the orbitals, and the columns represent
     * the positions.
     *
     * This matrix is initialized in the constructor of a subclass, and it
     * should remain constant after initialization is complete.
     */
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> orbitals;

    const boost::shared_ptr<const Lattice> lattice;
};

#endif
