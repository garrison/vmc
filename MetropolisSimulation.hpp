#ifndef _METROPOLIS_SIMULATION_HPP
#define _METROPOLIS_SIMULATION_HPP

#ifdef DEBUG
#include <iostream>
#endif

#include <boost/assert.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random.hpp>

#include "vmc-typedefs.hpp"

template <class Walk_T>
class MetropolisSimulation
{
public:
    MetropolisSimulation (const Walk_T &walk_, unsigned int lg_initialization_sweeps, const rng_seed_t &seed)
	: walk(walk_),
	  m_steps(0),
	  m_steps_accepted(0),
	  rng(seed),
	  uniform_distribution(rng_class(rng())) // see http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
	{
	    // do initialization sweeps
	    unsigned int initialization_sweeps = 1 << lg_initialization_sweeps;
	    for (unsigned int i = 0; i < initialization_sweeps; ++i)
		perform_single_step();
	}

    void iterate (unsigned int lg_sweeps)
	{
	    // fixme: steps between measurement, and implement measurement
	    unsigned int sweeps = 1 << lg_sweeps;
	    for (unsigned int i = 0; i < sweeps; ++i)
		perform_single_step();
	}

    int steps_completed (void) const
	{
	    return m_steps;
	}

    int steps_rejected (void) const
	{
	    return m_steps - m_steps_accepted;
	}

    int steps_accepted (void) const
	{
	    return m_steps_accepted;
	}

    const Walk_T & get_walk (void) const
	{
	    return walk;
	}

private:
    Walk_T walk;
    int m_steps, m_steps_accepted;

    rng_class rng;
    boost::uniform_01<rng_class> uniform_distribution;

    void perform_single_step (void)
	{
	    Walk_T proposed_step(walk);

	    probability_t probability_ratio = proposed_step.compute_probability_ratio_of_random_transition(rng);
	    BOOST_ASSERT(probability_ratio >= 0);
	    if (probability_ratio >= 1
		|| probability_ratio > uniform_distribution()) {
		// accept transition
		proposed_step.accept_transition();
		walk = proposed_step;
		++m_steps_accepted;
#ifdef DEBUG
		std::cerr << "A" << std::endl;
#endif
	    } else {
#ifdef DEBUG
		std::cerr << "-" << std::endl;
#endif
	    }
	    ++m_steps;
	}

};

#endif
