#ifndef _VMC_ORBITAL_DEFINITIONS_HPP
#define _VMC_ORBITAL_DEFINITIONS_HPP

#include <memory>
#include <cassert>

// we always want to include vmc-typedefs.hpp before including Eigen
#include "vmc-typedefs.hpp"
#include <Eigen/Dense>

class Lattice;

/**
 * Class for defining the orbitals in a determinantal wavefunction
 */
template <typename AmplitudeType>
class OrbitalDefinitions
{
public:
    OrbitalDefinitions (const Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> &orbitals_, const std::shared_ptr<const Lattice> &lattice_)
        : orbitals(orbitals_),
          lattice(lattice_)
        {
            assert(orbitals.cols() == lattice_->total_sites());
        }

#if 0
    virtual ~OrbitalDefinitions (void)
        {
        }
#endif

    /**
     * Returns the set of orbitals at the given position
     */
    typename Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic>::ConstColXpr at_position (unsigned int site_index) const
        {
            assert(site_index < orbitals.cols());
            return orbitals.col(site_index);
        }

    /**
     * Returns the matrix representing the orbitals
     *
     * The rows of the matrix represent the orbitals, and the columns represent
     * each possible position.
     */
    const Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> & get_orbitals (void) const
        {
            return orbitals;
        }

    unsigned int get_N_filled (void) const
        {
            return orbitals.rows();
        }

    unsigned int get_N_sites (void) const
        {
            return orbitals.cols();
        }

    const std::shared_ptr<const Lattice> & get_lattice_ptr (void) const
        {
            return lattice;
        }

    // CYTHON-LIMITATION: http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#c-left-values
    static inline void set_matrix_coeff (Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> &mat, unsigned int row, unsigned int col, AmplitudeType value)
        {
            mat(row, col) = value;
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
    Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> orbitals;

    const std::shared_ptr<const Lattice> lattice;
};

#endif
