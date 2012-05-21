#ifndef _GREEN_MEASUREMENT_HPP
#define _GREEN_MEASUREMENT_HPP

#include <Eigen/Core>
#include <boost/assert.hpp>
#include <boost/cast.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "NDLattice.hpp"
#include "PositionArguments.hpp"

/**
 * Green function measurement
 *
 * Assumes translational invariance
 *
 * @see StandardWalk
 */
template <std::size_t DIM>
class GreenMeasurement : public Measurement<StandardWalk>
{
public:
    GreenMeasurement (unsigned int steps_per_measurement, unsigned int species_)
        : Measurement<StandardWalk>(steps_per_measurement),
          species(species_),
          denominator(0)
        {
        }

    /**
     * Returns the Green function measurement so far for a given vector
     */
    amplitude_t get (unsigned int site_index, unsigned int basis_index=0) const
        {
            BOOST_ASSERT(site_index < get_N_sites());
            BOOST_ASSERT(basis_index < basis_indices());
            return green_accum(basis_index, site_index) / real_t(denominator);
        }

    /**
     * Number of sites per unit cell of the Bravais lattice
     */
    unsigned int basis_indices (void) const
        {
            return green_accum.rows();
        }

    /**
     * Number of sites on the lattice
     */
    unsigned int get_N_sites (void) const
        {
            return green_accum.cols();
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
            green_accum.setZero(basis_indices, total_sites);
            current_step_green_accum.resizeLike(green_accum);

            single_step_denominator = lattice->total_sites() / basis_indices;
        }

    /**
     * Calculate and tally a measurement
     */
    void measure_ (const StandardWalk &walk)
        {
            const WavefunctionAmplitude &wf = walk.get_wavefunction();
            const PositionArguments &r = wf.get_positions();
            const NDLattice<DIM> *lattice = boost::polymorphic_downcast<const NDLattice<DIM>*>(&wf.get_lattice());

            current_step_green_accum.setZero();

            // loop through all (particle, empty site) pairs
            for (unsigned int i = 0; i < r.get_N_filled(species); ++i) {
                const Particle particle(i, species);
                const typename NDLattice<DIM>::Site site_i(lattice->site_from_index(r[particle]));

                // amplitude of starting and ending on same site
                current_step_green_accum(site_i.basis_index, 0) += amplitude_t(1);

                // loop through all empty sites
                for (unsigned int j = 0; j < r.get_N_sites(); ++j) {
                    if (r.is_occupied(j, species))
                        continue;

                    boost::shared_ptr<WavefunctionAmplitude> wf_operated = wf.clone();
                    wf_operated->move_particle(particle, j);

                    typename NDLattice<DIM>::Site site_j(lattice->site_from_index(j));
                    phase_t phase = lattice->asm_subtract_site_vector(site_j, site_i.bravais_site());
                    // fixme: check logic of multiplying by phase
                    current_step_green_accum(site_i.basis_index, lattice->site_to_index(site_j)) += std::conj(wf_operated->psi() * phase / wf.psi());
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
            green_accum += current_step_green_accum;
            denominator += single_step_denominator;
        }

    const unsigned int species;

    // row is the basis, column is the site index
    Eigen::Array<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> green_accum, current_step_green_accum;

    unsigned int denominator, single_step_denominator;
};

#endif
