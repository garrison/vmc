#ifndef _DENSITY_DENSITY_MEASUREMENT_HPP
#define _DENSITY_DENSITY_MEASUREMENT_HPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"

/**
 * Density-density measurement
 *
 * Assumes translational invariance
 *
 * @see StandardWalk
 */
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
    real_t get (unsigned int site_index, unsigned int basis_index=0) const;

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
    void initialize_ (const StandardWalk &walk);

    /**
     * Calculate and tally a measurement
     */
    void measure_ (const StandardWalk &walk);

    /**
     * Tally again the most recent measurement
     */
    void repeat_measurement_ (const StandardWalk &walk);

    unsigned int species1, species2;

    // row is the basis, column is the site index
    Eigen::Array<unsigned int, Eigen::Dynamic, Eigen::Dynamic> density_accum, current_step_density_accum;

    unsigned int denominator, single_step_denominator;

    real_t density_squared;
};

#endif
