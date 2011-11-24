#ifndef _ORBITAL_DEFINITIONS_HPP
#define _ORBITAL_DEFINITIONS_HPP

#include <Eigen/Dense>
#include <boost/assert.hpp>

#include "vmc-typedefs.hpp"
#include "Lattice.hpp"

class OrbitalDefinitions
// abstract base class
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

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic>::ConstColXpr get (unsigned int site_index) const
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
    // this should remain constant after initialization is complete
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> orbitals;

    const boost::shared_ptr<const Lattice> lattice;
};

#endif
