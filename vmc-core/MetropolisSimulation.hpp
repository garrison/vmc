#ifndef _METROPOLIS_SIMULATION_HPP
#define _METROPOLIS_SIMULATION_HPP

#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
#include <iostream>
#endif

#include <list>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random.hpp>

#include "vmc-typedefs.hpp"
#include "Measurement.hpp"

/**
 * Metropolis Simulation
 *
 * This class performs the Metropolis Monte Carlo algorithm using some given
 * walk.  First it does some number of steps to (hopefully) reach
 * "equilibrium," and after that it does some more steps, taking a measurement
 * after each move.
 */
template <class Walk_T>
class MetropolisSimulation
{
public:
    /**
     * Constructor
     *
     * Note: the initialization sweeps are performed during object
     * construction.
     *
     * @param walk_ walk in its initial state
     *
     * @param measurements_ a list of measurements to perform
     *
     * @param initialization_sweeps number of steps to take before beginning to
     * take measurements
     *
     * @param seed random seed
     *
     * The Walk_T type must define three methods:
     *
     * * probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
     * * void accept_transition ();
     * * void reject_transition ();
     *
     * @see StandardWalk for an example walk
     */
    MetropolisSimulation (const Walk_T &walk_, const std::list<boost::shared_ptr<Measurement<Walk_T> > > &measurements_,
                          unsigned int initialization_sweeps, const rng_seed_t &seed)
        : walk(walk_),
          measurements(measurements_),
          m_steps(0),
          m_steps_accepted(0),
          m_steps_fully_rejected(0),
          measurement_not_yet_updated(true),
          rng(seed),
          uniform_distribution(rng_class(rng())) // see http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
        {
            perform_initialization(initialization_sweeps);
        }

    /**
     * Convenience constructor, if only one measurement is being performed
     *
     * The initialization sweeps are performed during object construction.
     *
     * @param walk_ walk in its initial state
     *
     * @param measurement_ the measurement to perform
     *
     * @param initialization_sweeps number of steps to take before beginning to
     * take measurements
     *
     * @param seed random seed
     */
    MetropolisSimulation (const Walk_T &walk_, const boost::shared_ptr<Measurement<Walk_T> > &measurement_,
                          unsigned int initialization_sweeps, const rng_seed_t &seed)
        : walk(walk_),
          m_steps(0),
          m_steps_accepted(0),
          m_steps_fully_rejected(0),
          measurement_not_yet_updated(true),
          rng(seed),
          uniform_distribution(rng_class(rng())) // see http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
        {
            measurements.push_back(measurement_);
            perform_initialization(initialization_sweeps);
        }

    /**
     * Perform some number of steps on the system, taking a measurement each
     * time
     */
    void iterate (unsigned int sweeps)
        {
            for (unsigned int i = 0; i < sweeps; ++i) {
                // perform a single step
                const bool accepted = perform_single_step();

                // perform the measurement(s)
                if (accepted || measurement_not_yet_updated) {
                    for (typename std::list<boost::shared_ptr<Measurement<Walk_T> > >::iterator m = measurements.begin(); m != measurements.end(); ++m)
                        (*m)->step_advanced(walk);
                    measurement_not_yet_updated = false;
                } else {
                    for (typename std::list<boost::shared_ptr<Measurement<Walk_T> > >::iterator m = measurements.begin(); m != measurements.end(); ++m)
                        (*m)->step_repeated(walk);
                }
            }
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

    /**
     * Returns the walk object
     */
    const Walk_T & get_walk (void) const
        {
            return walk;
        }

private:
    Walk_T walk;
    std::list<boost::shared_ptr<Measurement<Walk_T> > > measurements;
    unsigned int m_steps, m_steps_accepted, m_steps_fully_rejected;
    bool measurement_not_yet_updated;

    rng_class rng;
    boost::uniform_01<rng_class> uniform_distribution;

    bool perform_single_step (void)
        {
            probability_t probability_ratio = walk.compute_probability_ratio_of_random_transition(rng);
            ++m_steps;
#if 1
            if (!(probability_ratio >= 0))
                std::cerr << "invalid probability ratio: " << probability_ratio << std::endl;
#endif
            if (probability_ratio >= 1
                || (probability_ratio > 0 && probability_ratio > uniform_distribution())) {
                // accept transition
#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
                std::cerr << "A" << std::endl;
#endif
                walk.accept_transition();
                ++m_steps_accepted;
                return true;
            } else {
                // reject transition
#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
                std::cerr << "-" << std::endl;
#endif
                walk.reject_transition();
                if (probability_ratio == 0)
                    ++m_steps_fully_rejected;
                return false;
            }
        }

    void perform_initialization (unsigned int initialization_sweeps)
        {
            // do initialization sweeps
            for (unsigned int i = 0; i < initialization_sweeps; ++i)
                perform_single_step();
            // initialize the measurements
            for (typename std::list<boost::shared_ptr<Measurement<Walk_T> > >::iterator m = measurements.begin(); m != measurements.end(); ++m)
                (*m)->initialize(walk);
        }

};

#endif
