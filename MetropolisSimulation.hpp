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

template <class Walk_T>
class MetropolisSimulation
{
public:
    MetropolisSimulation (const Walk_T &walk_, const std::list<boost::shared_ptr<Measurement<Walk_T> > > &measurements_,
                          unsigned int lg_initialization_sweeps, const rng_seed_t &seed)
        : walk(walk_),
          measurements(measurements_),
          m_steps(0),
          m_steps_accepted(0),
          m_steps_fully_rejected(0),
          not_yet_measured(true),
          rng(seed),
          uniform_distribution(rng_class(rng())) // see http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
        {
            perform_initialization(lg_initialization_sweeps);
        }

    MetropolisSimulation (const Walk_T &walk_, const boost::shared_ptr<Measurement<Walk_T> > &measurement_,
                          unsigned int lg_initialization_sweeps, const rng_seed_t &seed)
        : walk(walk_),
          m_steps(0),
          m_steps_accepted(0),
          m_steps_fully_rejected(0),
          not_yet_measured(true),
          rng(seed),
          uniform_distribution(rng_class(rng())) // see http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
        {
            measurements.push_back(measurement_);
            perform_initialization(lg_initialization_sweeps);
        }

    void iterate (unsigned int lg_sweeps)
        {
            unsigned int sweeps = 1 << lg_sweeps;
            for (unsigned int i = 0; i < sweeps; ++i) {
                // perform a single step
                const bool accepted = perform_single_step();

                // perform the measurement(s)
                if (accepted || not_yet_measured) {
                    for (typename std::list<boost::shared_ptr<Measurement<Walk_T> > >::iterator m = measurements.begin(); m != measurements.end(); ++m)
                        (*m)->measure(walk);
                    not_yet_measured = false;
                } else {
                    for (typename std::list<boost::shared_ptr<Measurement<Walk_T> > >::iterator m = measurements.begin(); m != measurements.end(); ++m)
                        (*m)->repeat_measurement(walk);
                }
            }
        }

    unsigned int steps_completed (void) const
        {
            return m_steps;
        }

    unsigned int steps_rejected (void) const
        {
            return m_steps - m_steps_accepted;
        }

    unsigned int steps_fully_rejected (void) const
        {
            return m_steps_fully_rejected;
        }

    unsigned int steps_accepted (void) const
        {
            return m_steps_accepted;
        }

    const Walk_T & get_walk (void) const
        {
            return walk;
        }

private:
    Walk_T walk;
    std::list<boost::shared_ptr<Measurement<Walk_T> > > measurements;
    unsigned int m_steps, m_steps_accepted, m_steps_fully_rejected;
    bool not_yet_measured;

    rng_class rng;
    boost::uniform_01<rng_class> uniform_distribution;

    bool perform_single_step (void)
        {
            Walk_T proposed_step(walk);

            probability_t probability_ratio = proposed_step.compute_probability_ratio_of_random_transition(rng);
            ++m_steps;
#if 1
            if (!(probability_ratio >= 0))
                std::cerr << "invalid probability ratio: " << probability_ratio << std::endl;
#endif
            if (probability_ratio >= 1
                || (probability_ratio > 0 && probability_ratio > uniform_distribution())) {
                // accept transition
                proposed_step.accept_transition();
                walk = proposed_step;
                ++m_steps_accepted;
#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
                std::cerr << "A" << std::endl;
#endif
                return true;
            } else {
                if (probability_ratio == 0)
                    ++m_steps_fully_rejected;
#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
                std::cerr << "-" << std::endl;
#endif
                return false;
            }
        }

    void perform_initialization (unsigned int lg_initialization_sweeps)
        {
            // do initialization sweeps
            unsigned int initialization_sweeps = 1 << lg_initialization_sweeps;
            for (unsigned int i = 0; i < initialization_sweeps; ++i)
                perform_single_step();
            // initialize the measurements
            for (typename std::list<boost::shared_ptr<Measurement<Walk_T> > >::iterator m = measurements.begin(); m != measurements.end(); ++m)
                (*m)->initialize(walk);
        }

};

#endif
