#ifndef _METROPOLIS_SIMULATION_HPP
#define _METROPOLIS_SIMULATION_HPP

#include <iostream>
#include <boost/assert.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random.hpp>

#include "vmc-typedefs.hpp"

template <class WaveFunction_T>
class MetropolisSimulation
{
public:
    MetropolisSimulation (const WaveFunction_T &wf, unsigned int lg_initialization_sweeps)
	: wavefunction(wf),
	  arguments(*WaveFunction_T::ArgumentWalk::random_initial_state(wf)),
	  m_steps(0),
	  m_steps_accepted(0),
	  rng(0), // fixme: seed
	  uniform_distribution(rng)
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

private:
    const WaveFunction_T wavefunction;
    typename WaveFunction_T::ArgumentWalk arguments;
    int m_steps, m_steps_accepted;

    // see http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
    boost::mt19937 rng;
    boost::uniform_01<boost::mt19937> uniform_distribution;

    void perform_single_step (void)
	{
	    typename WaveFunction_T::ArgumentWalk new_arguments(arguments);

	    probability_t probability_ratio = new_arguments.compute_probability_ratio_of_random_transition();
	    BOOST_ASSERT(probability_ratio >= 0);
	    if (probability_ratio >= 1
		|| probability_ratio > uniform_distribution()) {
		// accept transition
		new_arguments.finalize_transition();
		arguments = new_arguments;
		++m_steps_accepted;
		std::cerr << "A" << std::endl;
	    } else {
		std::cerr << "-" << std::endl;
	    }
	    ++m_steps;
	}

};

#endif
