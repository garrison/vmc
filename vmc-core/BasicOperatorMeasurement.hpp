#ifndef _VMC_BASIC_OPERATOR_MEASUREMENT_HPP
#define _VMC_BASIC_OPERATOR_MEASUREMENT_HPP

#include <memory>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "BasicOperator.hpp"
#include "BasicOperatorEvaluator.hpp"
#include "BlockedEstimate.hpp"

/**
 * Basic operator measurement
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
class BasicOperatorMeasurement : public Measurement<StandardWalk<_AmplitudeType> >
{
public:
    typedef _AmplitudeType AmplitudeType;
    typedef AmplitudeType PhaseType;
    typedef Walk<typename StandardWalk<AmplitudeType>::ProbabilityType> BaseWalkType;
    typedef StandardWalk<AmplitudeType> WalkType;

    BasicOperatorMeasurement (unsigned int steps_per_measurement,
                              const BasicOperator &operator_,
                              const BoundaryConditions<PhaseType> &bcs_)
        : Measurement<WalkType>(steps_per_measurement),
          evaluator(operator_),
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
     * Calculate and tally a measurement
     */
    virtual void measure_ (const WalkType &walk) override
        {
            const typename Wavefunction<AmplitudeType>::Amplitude &wfa = walk.get_wavefunctionamplitude();
            most_recent_value = evaluator.evaluate(wfa, bcs);
            estimate.add_value(most_recent_value);
        }

    /**
     * Tally again the most recent measurement
     */
    virtual void repeat_measurement_ (const WalkType &walk) override
        {
            (void) walk;
            estimate.add_value(most_recent_value);
        }

    virtual bool is_valid_walk_ (const WalkType &walk) override
        {
            return BasicOperator::is_valid(evaluator.m_operator.hopv,
                                           walk.get_wavefunctionamplitude().get_lattice(),
                                           walk.get_wavefunctionamplitude().get_positions().get_N_species());
        }

    const BasicOperatorEvaluator evaluator;
    const BoundaryConditions<PhaseType> bcs;

    BlockedEstimate<AmplitudeType> estimate;
    AmplitudeType most_recent_value;
};

#endif
