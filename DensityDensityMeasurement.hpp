#ifndef _DENSITY_DENSITY_MEASUREMENT_HPP
#define _DENSITY_DENSITY_MEASUREMENT_HPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "Lattice.hpp"
#include "PositionArguments.hpp"

template <class Lattice_T>
class DensityDensityMeasurement : public Measurement<StandardWalk>
{
public:
    real_t get (unsigned int site_index, unsigned int basis_index=0) const
	{
	    BOOST_ASSERT(site_index < density_accum.cols());
	    BOOST_ASSERT(basis_index < density_accum.rows());
	    unsigned int num = density_accum(basis_index, site_index);
	    return real_t(num) / denominator(basis_index);
	}

private:
    void initialize_ (const StandardWalk &walk)
	{
	    const unsigned int total_sites = walk.get_wavefunction().get_lattice().total_sites();
	    BOOST_ASSERT(total_sites > 0);
	    const Lattice_T *lattice = dynamic_cast<const Lattice_T*>(&walk.get_wavefunction().get_lattice());
	    BOOST_ASSERT(lattice != 0);

	    const unsigned int basis_indices = lattice->basis_indices;
	    density_accum.setZero(basis_indices, total_sites);
	    denominator.setZero(basis_indices);
	}

    void measure_ (const StandardWalk &walk)
	{
	    const PositionArguments &r = walk.get_wavefunction().get_positions();
	    const Lattice_T *lattice = dynamic_cast<const Lattice_T*>(&walk.get_wavefunction().get_lattice());
	    BOOST_ASSERT(lattice != 0);

	    const unsigned int total_sites = walk.get_wavefunction().get_lattice().total_sites();
	    const unsigned int basis_indices = lattice->basis_indices;
	    current_density_accum.setZero(basis_indices, total_sites);
	    current_denominator.setZero(basis_indices);

	    // loop through all pairs of particles
	    for (unsigned int i = 0; i < r.get_N_filled(); ++i) {
		typename Lattice_T::Site site_i(lattice->site_from_index(r[i]));
		unsigned int i_basis = lattice->basis_index(site_i);
		site_i = lattice->move_to_basis_index(site_i, 0);
		for (unsigned int j = 0; j < r.get_N_filled(); ++j) {
		    typename Lattice_T::Site site_j(lattice->site_from_index(r[j]));
		    lattice->asm_subtract_site_vector(site_j, site_i);
		    ++density_accum(i_basis, lattice->site_to_index(site_j));
		}
		++denominator(i_basis);
	    }

	    repeat_measurement_(walk);
	}

    void repeat_measurement_ (const StandardWalk &walk)
	{
	    (void) walk;
	    density_accum += current_density_accum;
	    denominator += current_denominator;
	}

    // row is the basis, column is the site index
    Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic> density_accum, current_density_accum;

    // index refers to the basis
    Eigen::Array<unsigned int, Eigen::Dynamic, 1> denominator, current_denominator;
};

#endif
