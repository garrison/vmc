#ifndef _VMC_OPERATOR_MEASUREMENT_HPP
#define _VMC_OPERATOR_MEASUREMENT_HPP

#include <memory>

#include <boost/assert.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "BasicOperator.hpp"
#include "BlockedEstimate.hpp"

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
template <typename _AmplitudeType>
class OperatorMeasurement : public Measurement<StandardWalk<_AmplitudeType> >
{
public:
    typedef _AmplitudeType AmplitudeType;
    typedef AmplitudeType PhaseType;
    typedef Walk<typename StandardWalk<AmplitudeType>::ProbabilityType> BaseWalkType;
    typedef StandardWalk<AmplitudeType> WalkType;

    OperatorMeasurement (unsigned int steps_per_measurement,
                         const BasicOperator &operator_,
                         const BoundaryConditions<PhaseType> &bcs_)
        : Measurement<WalkType>(steps_per_measurement),
          m_operator(operator_),
          bcs(bcs_)
        {
        }

    /**
     * Returns the operator measurement estimate for a given vector
     */
    const BlockedEstimate<AmplitudeType> & get_estimate (void) const
        {
            return estimate;
        }

private:
    /**
     * Prepare the object for taking measurements
     */
    virtual void initialize_ (const WalkType &walk) override;

    /**
     * Calculate and tally a measurement
     */
    virtual void measure_ (const WalkType &walk) override;

    /**
     * Tally again the most recent measurement
     */
    virtual void repeat_measurement_ (const WalkType &walk) override;

    virtual bool is_valid_walk_ (const WalkType &walk) override
        {
            return BasicOperator::is_valid(m_operator.hopv,
                                           walk.get_wavefunctionamplitude().get_lattice(),
                                           walk.get_wavefunctionamplitude().get_positions().get_N_species());
        }

    bool is_sum_over_sites (void) const
        {
            return bcs.size() != 0;
        }

    const BasicOperator m_operator;
    const BoundaryConditions<PhaseType> bcs;

    BlockedEstimate<AmplitudeType> estimate;
    AmplitudeType most_recent_value;
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class OperatorMeasurement<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif
