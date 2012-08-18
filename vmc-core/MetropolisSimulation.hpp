#ifndef _METROPOLIS_SIMULATION_HPP
#define _METROPOLIS_SIMULATION_HPP

#include <memory>
#include <list>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include "vmc-typedefs.hpp"
#include "RandomNumberGenerator.hpp"

class BaseMeasurement;
class Walk;

/**
 * Metropolis Simulation
 *
 * This class performs the Metropolis Monte Carlo algorithm using some given
 * walk.  First it does some number of steps to (hopefully) reach
 * "equilibrium," and after that it does some more steps, taking a measurement
 * after each move.
 */
class MetropolisSimulation : boost::noncopyable
{
public:
    /**
     * Constructor
     *
     * Note: the initialization sweeps are performed during object
     * construction.
     *
     * @param walk_ walk in its initial state; since it's an auto_ptr, the new
     * MetropolisSimulation object will be passed ownership of the Walk
     *
     * @param measurements_ a list of measurements to perform
     *
     * @param initialization_sweeps number of steps to take before beginning to
     * take measurements
     *
     * @param rng random number generator; since it's an auto_ptr, ownership
     * will be passed to the new MetropolisSimulation
     *
     * @see Walk
     */
    MetropolisSimulation (std::auto_ptr<Walk> &walk_, const std::list<boost::shared_ptr<BaseMeasurement> > &measurements_,
                          unsigned int initialization_sweeps, std::auto_ptr<RandomNumberGenerator> &rng_);

    /**
     * Perform some number of steps on the system, taking a measurement each
     * time
     */
    void iterate (unsigned int sweeps);

    /**
     * Returns the walk object
     */
    const Walk * get_walk_ptr (void) const
        {
            return walk.get();
        }

    /**
     * Returns the number of steps completed so far
     */
    unsigned int steps_completed (void) const
        {
            return m_steps;
        }

    /**
     * Returns the number of rejected moves so far
     */
    unsigned int steps_rejected (void) const
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
    unsigned int steps_fully_rejected (void) const
        {
            return m_steps_fully_rejected;
        }

    /**
     * Returns the number of accepted moves so far
     */
    unsigned int steps_accepted (void) const
        {
            return m_steps_accepted;
        }

    const std::list<boost::shared_ptr<BaseMeasurement> > & get_measurements (void) const
        {
            // fixme: in theory, somebody could do something evil and directly
            // modify the measurements returned here, since they aren't
            // const-protected.  but, they could also keep references to the
            // measurements since they were passed in the constructor, so it
            // makes no difference at the moment.
            return measurements;
        }

protected:
    unsigned int m_steps, m_steps_accepted, m_steps_fully_rejected;

private:
    std::auto_ptr<Walk> walk;
    std::list<boost::shared_ptr<BaseMeasurement> > measurements;
    bool measurement_not_yet_updated;

    std::auto_ptr<RandomNumberGenerator> rng;

    bool perform_single_step (void);

    void perform_initialization (unsigned int initialization_sweeps);
};

#endif
