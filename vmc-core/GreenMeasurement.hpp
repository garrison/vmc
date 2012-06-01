#ifndef _GREEN_MEASUREMENT_HPP
#define _GREEN_MEASUREMENT_HPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"

/**
 * Green function measurement
 *
 * Assumes translational invariance
 *
 * @see StandardWalk
 */
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
    void initialize_ (const StandardWalk &walk);

    /**
     * Calculate and tally a measurement
     */
    void measure_ (const StandardWalk &walk);

    /**
     * Tally again the most recent measurement
     */
    void repeat_measurement_ (const StandardWalk &walk);

    const unsigned int species;

    // row is the basis, column is the site index
    Eigen::Array<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> green_accum, current_step_green_accum;

    unsigned int denominator, single_step_denominator;
};

#endif
