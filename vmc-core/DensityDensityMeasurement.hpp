#ifndef _DENSITY_DENSITY_MEASUREMENT_HPP
#define _DENSITY_DENSITY_MEASUREMENT_HPP

#include <Eigen/Core>
#include <boost/assert.hpp>
#include <boost/cast.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "NDLattice.hpp"
#include "PositionArguments.hpp"

/**
 * Density-density measurement
 *
 * Assumes translational invariance
 *
 * @see StandardWalk
 */
template <std::size_t DIM>
class DensityDensityMeasurement : public Measurement<StandardWalk>
{
public:
    DensityDensityMeasurement (unsigned int steps_per_measurement, unsigned int species1_, unsigned int species2_)
        : Measurement<StandardWalk>(steps_per_measurement),
          species1(species1_),
          species2(species2_),
          denominator(0)
        {
        }

    /**
     * Returns the density-density measurement so far for a given vector
     */
    real_t get (unsigned int site_index, unsigned int basis_index=0) const
        {
            BOOST_ASSERT(site_index < density_accum.cols());
            BOOST_ASSERT(basis_index < density_accum.rows());
            unsigned int num = density_accum(basis_index, site_index);
            return real_t(num) / denominator - density_squared;
        }

    /**
     * Number of sites per unit cell of the Bravais lattice
     */
    unsigned int basis_indices (void) const
        {
            return density_accum.rows();
        }

    /**
     * Number of sites on the lattice
     */
    unsigned int get_N_sites (void) const
        {
            return density_accum.cols();
        }

private:
    /**
     * Prepare the object for taking measurements
     */
    void initialize_ (const StandardWalk &walk)
        {
            const unsigned int total_sites = walk.get_wavefunction().get_lattice().total_sites();
            BOOST_ASSERT(total_sites > 0);
            const NDLattice<DIM> *lattice = boost::polymorphic_downcast<const NDLattice<DIM>*>(&walk.get_wavefunction().get_lattice());

            const unsigned int basis_indices = lattice->basis_indices;
            density_accum.setZero(basis_indices, total_sites);
            current_step_density_accum.resizeLike(density_accum);

            single_step_denominator = lattice->total_sites() / basis_indices;

            const PositionArguments &r = walk.get_wavefunction().get_positions();

            BOOST_ASSERT(species1 < r.get_N_species());
            BOOST_ASSERT(species2 < r.get_N_species());

            density_squared = real_t(r.get_N_filled(species1) * r.get_N_filled(species2)) / (r.get_N_sites() * r.get_N_sites());
        }

    /**
     * Calculate and tally a measurement
     */
    void measure_ (const StandardWalk &walk)
        {
            const PositionArguments &r = walk.get_wavefunction().get_positions();
            const NDLattice<DIM> *lattice = boost::polymorphic_downcast<const NDLattice<DIM>*>(&walk.get_wavefunction().get_lattice());

            current_step_density_accum.setZero();

            // loop through all pairs of particles
            for (unsigned int i = 0; i < r.get_N_filled(species2); ++i) {
                const typename NDLattice<DIM>::Site site_i(lattice->site_from_index(r[Particle(i, species2)]));
                for (unsigned int j = 0; j < r.get_N_filled(species1); ++j) {
                    typename NDLattice<DIM>::Site site_j(lattice->site_from_index(r[Particle(j, species1)]));
                    lattice->asm_subtract_site_vector(site_j, site_i.bravais_site());
                    ++current_step_density_accum(site_i.basis_index, lattice->site_to_index(site_j));
                }
            }

            repeat_measurement_(walk);
        }

    /**
     * Tally again the most recent measurement
     */
    void repeat_measurement_ (const StandardWalk &walk)
        {
            (void) walk;
            density_accum += current_step_density_accum;
            denominator += single_step_denominator;
        }

    unsigned int species1, species2;

    // row is the basis, column is the site index
    Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic> density_accum, current_step_density_accum;

    unsigned int denominator, single_step_denominator;

    real_t density_squared;
};

#endif
