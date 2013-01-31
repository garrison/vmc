#ifndef _VMC_OPERATOR_MEASUREMENT_HPP
#define _VMC_OPERATOR_MEASUREMENT_HPP

#include <memory>

#include <boost/assert.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "ParticleOperator.hpp"
#include "BinnedEstimate.hpp"

/**
 * Operator measurement
 *
 * Measures the given operator on the lattice.
 *
 * If boundary conditions are given, this measurement will loop over all sites,
 * summing the expectation value of the operator as translated across the
 * lattice (wrapping around any boundary conditions that are not open).  If
 * `boundary_conditions` is the empty array, no sum is performed.
 *
 * @see StandardWalk
 */
class OperatorMeasurement : public Measurement<StandardWalk>
{
public:
    OperatorMeasurement (unsigned int steps_per_measurement,
                         const ParticleOperator &operator_,
                         const BoundaryConditions &bcs_)
        : Measurement<StandardWalk>(steps_per_measurement),
          m_operator(operator_),
          bcs(bcs_)
        {
        }

    /**
     * Returns the operator measurement estimate for a given vector
     */
    const BinnedEstimate<amplitude_t> & get_estimate (void) const
        {
            return estimate;
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

    bool is_valid_walk_ (const StandardWalk &walk)
        {
            return ParticleOperator::is_valid(m_operator.hopv,
                                              walk.get_wavefunctionamplitude().get_lattice(),
                                              walk.get_wavefunctionamplitude().get_positions().get_N_species());
        }

    bool is_sum_over_sites (void) const
        {
            return bcs.size() != 0;
        }

    const ParticleOperator m_operator;
    const BoundaryConditions bcs;

    BinnedEstimate<amplitude_t> estimate;
    amplitude_t most_recent_value;
};

#endif
