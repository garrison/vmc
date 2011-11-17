#ifndef _DENSITY_DENSITY_MEASUREMENT_HPP
#define _DENSITY_DENSITY_MEASUREMENT_HPP

#include <vector>

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
	    BOOST_ASSERT(site_index < density_accum[0].size());
	    BOOST_ASSERT(basis_index < density_accum.size());
	    unsigned int num = density_accum[basis_index][site_index];
	    return real_t(num) / denominator[basis_index];
	}

private:
    void initialize_ (const StandardWalk &walk)
	{
	    const unsigned int total_sites = walk.get_wavefunction().get_lattice().total_sites();
	    BOOST_ASSERT(total_sites > 0);
	    const Lattice_T *lattice = dynamic_cast<const Lattice_T*>(&walk.get_wavefunction().get_lattice());
	    BOOST_ASSERT(lattice != 0);

	    const unsigned int basis_indices = lattice->basis_indices;
	    density_accum.resize(basis_indices);
	    denominator.resize(basis_indices);
	    for (unsigned int i = 0; i < basis_indices; ++i)
		density_accum[i].resize(total_sites);
	}

    void measure_ (const StandardWalk &walk)
	{
	    const PositionArguments &r = walk.get_wavefunction().get_positions();
	    const Lattice_T *lattice = dynamic_cast<const Lattice_T*>(&walk.get_wavefunction().get_lattice());
	    BOOST_ASSERT(lattice != 0);

	    // loop through all pairs of particles
	    for (unsigned int i = 0; i < r.get_N_filled(); ++i) {
		typename Lattice_T::Site site_i(lattice->site_from_index(r[i]));
		unsigned int i_basis = lattice->basis_index(site_i);
		site_i = lattice->move_to_basis_index(site_i, 0);
		for (unsigned int j = 0; j < r.get_N_filled(); ++j) {
		    typename Lattice_T::Site site_j(lattice->site_from_index(r[j]));
		    lattice->asm_subtract_site_vector(site_j, site_i);
		    ++density_accum[i_basis][lattice->site_to_index(site_j)];
		}
		++denominator[i_basis];
	    }
	}

    // first index is the basis, second is the site index
    std::vector<std::vector<unsigned int> > density_accum;

    // index refers to the basis
    std::vector<unsigned int> denominator;
};

#endif
