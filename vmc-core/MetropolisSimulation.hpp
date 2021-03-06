#ifndef _VMC_METROPOLIS_SIMULATION_HPP
#define _VMC_METROPOLIS_SIMULATION_HPP

#include <memory>
#include <list>
#include <cassert>

#include <boost/noncopyable.hpp>

#include "vmc-typedefs.hpp"
#include "RandomNumberGenerator.hpp"

template <typename ProbabilityType>
class BaseMeasurement;

template <typename ProbabilityType>
class Walk;

/**
 * Metropolis Simulation
 *
 * This class performs the Metropolis Monte Carlo algorithm using some given
 * walk.  First it does some number of steps to (hopefully) reach
 * "equilibrium," and after that it does some more steps, taking a measurement
 * after each move.
 */
template <typename _ProbabilityType, typename _CounterType=unsigned long long>
class MetropolisSimulation : boost::noncopyable
{
public:
    typedef _ProbabilityType ProbabilityType;
    typedef _CounterType CounterType;
    typedef Walk<ProbabilityType> WalkType;
    typedef BaseMeasurement<ProbabilityType> BaseMeasurementType;

    /**
     * Constructor
     *
     * Note: the initialization sweeps are performed during object
     * construction.
     *
     * @param walk_ walk in its initial state; since it's a unique_ptr, the new
     * MetropolisSimulation object will be passed ownership of the Walk
     *
     * @param measurements_ a list of measurements to perform
     *
     * @param initialization_sweeps number of steps to take before beginning to
     * take measurements
     *
     * @param rng random number generator; since it's a unique_ptr, ownership
     * will be passed to the new MetropolisSimulation
     *
     * @see Walk
     */
    MetropolisSimulation (std::unique_ptr<WalkType> walk_, const std::list<std::shared_ptr<BaseMeasurementType> > &measurements_,
                          CounterType initialization_sweeps, std::unique_ptr<RandomNumberGenerator> rng_);

    /**
     * Perform some number of steps on the system, taking a measurement each
     * time
     */
    void iterate (CounterType sweeps);

    /**
     * Returns the walk object
     */
    const WalkType & get_walk (void) const
        {
            return *walk;
        }

    /**
     * Returns the number of steps completed so far
     */
    CounterType steps_completed (void) const
        {
            return m_steps;
        }

    /**
     * Returns the number of rejected moves so far
     */
    CounterType steps_rejected (void) const
        {
            return m_steps - m_steps_accepted;
        }

    /**
     * Returns the number of "fully rejected" moves, i.e. the number of
     * attempted moves that return a probability of precisely zero.
     *
     * This may be useful information to us if we are trying to improve the
     * statistics of the simulation.  If a large proportion of the moves
     * attempted are fully rejected, then maybe we should be choosing our moves
     * more wisely.
     *
     * One example of where this frequent rejection might occur is in a
     * wavefunction that disallows double-rung occupancy.  In such a case, it
     * would be smart for us not to even attempt moves that would lead to
     * states with zero amplitude.
     */
    CounterType steps_fully_rejected (void) const
        {
            return m_steps_fully_rejected;
        }

    /**
     * Returns the number of accepted moves so far
     */
    CounterType steps_accepted (void) const
        {
            return m_steps_accepted;
        }

    const std::list<std::shared_ptr<BaseMeasurementType> > & get_measurements (void) const
        {
            // fixme: in theory, somebody could do something evil and directly
            // modify the measurements returned here, since they aren't
            // const-protected.  but, they could also keep references to the
            // measurements since they were passed in the constructor, so it
            // makes no difference at the moment.
            return measurements;
        }

    /**
     * Performs any possible checks to see whether the simulation currently has
     * significant numerical error, and raises an exception if it does.
     */
    void check_for_numerical_error (void) const
        {
            walk->check_for_numerical_error();
        }

protected:
    CounterType m_steps, m_steps_accepted, m_steps_fully_rejected;

private:
    std::unique_ptr<WalkType> walk;
    std::list<std::shared_ptr<BaseMeasurementType> > measurements;
    bool measurement_not_yet_updated;

    std::unique_ptr<RandomNumberGenerator> rng;

    bool perform_single_step (void);

    void perform_initialization (CounterType initialization_sweeps);
};

#define VMC_SUPPORTED_REAL_TYPE(type) extern template class MetropolisSimulation<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_REAL_TYPE

#endif
