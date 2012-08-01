#ifndef _VMC_OPERATOR_MEASUREMENT_HPP
#define _VMC_OPERATOR_MEASUREMENT_HPP

#include <memory>

#include <boost/assert.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "ParticleOperator.hpp"
#include "BinnedEstimate.hpp"

template <class T>
static inline std::auto_ptr<T> ptr_to_auto_ptr (T *p)
{
    return std::auto_ptr<T>(p ? new T(*p) : 0);
}

/**
 * Operator measurement
 *
 * Measures the given operator on the lattice.
 *
 * If `sum` is true, it will loop over all sites, summing the expectation value
 * of the operator as translated across the lattice.  Further, if boundary
 * conditions are given, the operator sum will also count operators in the sum
 * that go off the lattice, there assuming that the lattice wraps with the
 * boundary conditions given.
 *
 * @see StandardWalk
 */
class OperatorMeasurement : public Measurement<StandardWalk>
{
public:
    OperatorMeasurement (unsigned int steps_per_measurement,
                         const ParticleOperator &operator_,
                         bool sum_=false,
                         const BoundaryConditions *bcs_=0)
        : Measurement<StandardWalk>(steps_per_measurement),
          m_operator(operator_),
          sum(sum_),
          bcs(ptr_to_auto_ptr(bcs_))
        {
            // boundary conditions should only be given if (sum == true)
            BOOST_ASSERT(!bcs_ || sum_);
        }

    /**
     * Returns the operator measurement so far for a given vector
     */
    amplitude_t get (void) const
        {
            return estimate.get_result();
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

    const ParticleOperator m_operator;
    bool sum;
    const std::auto_ptr<const BoundaryConditions> bcs;

    BinnedEstimate<amplitude_t> estimate;
    amplitude_t most_recent_value;
};

#endif
