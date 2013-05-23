#include <iostream>

// these three things are included for the exception handling
#include <stdexcept>
#include <string>
#include <boost/lexical_cast.hpp>

#include "Measurement.hpp"
#include "Walk.hpp"
#include "MetropolisSimulation.hpp"

template <typename ProbabilityType>
class invalid_probability_error : public std::range_error
{
public:
    invalid_probability_error (const ProbabilityType &invalid_probability_);

    static inline std::string construct_what_string (const ProbabilityType &invalid_probability)
        {
            return std::string("Invalid probability ratio: ") + boost::lexical_cast<std::string>(invalid_probability);
        }

private:
    const ProbabilityType invalid_probability;
};

template <typename ProbabilityType>
invalid_probability_error<ProbabilityType>::invalid_probability_error (const ProbabilityType &invalid_probability_)
    : std::range_error(construct_what_string(invalid_probability_)),
      invalid_probability(invalid_probability_)
{
}

template <typename ProbabilityType>
MetropolisSimulation<ProbabilityType>::MetropolisSimulation (std::unique_ptr<Walk<ProbabilityType> > &walk_, const std::list<boost::shared_ptr<BaseMeasurement<ProbabilityType> > > &measurements_,
                                            unsigned int initialization_sweeps, std::unique_ptr<RandomNumberGenerator> &rng_)
    : m_steps(0),
      m_steps_accepted(0),
      m_steps_fully_rejected(0),
      walk(std::move(walk_)),
      measurements(measurements_),
      measurement_not_yet_updated(true),
      rng(std::move(rng_))
{
    BOOST_ASSERT(rng.get() != 0);
#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    for (auto m = measurements.begin(); m != measurements.end(); ++m)
        BOOST_ASSERT((*m)->is_valid_walk(*walk));
#endif

    perform_initialization(initialization_sweeps);
}

template <typename ProbabilityType>
void MetropolisSimulation<ProbabilityType>::iterate (unsigned int sweeps)
{
    for (unsigned int i = 0; i < sweeps; ++i) {
        // perform a single step
        const bool accepted = perform_single_step();

        // perform the measurement(s)
        if (accepted || measurement_not_yet_updated) {
            for (auto m = measurements.begin(); m != measurements.end(); ++m)
                (*m)->step_advanced(*walk);
            measurement_not_yet_updated = false;
        } else {
            for (auto m = measurements.begin(); m != measurements.end(); ++m)
                (*m)->step_repeated(*walk);
        }
    }
}

template <typename ProbabilityType>
bool MetropolisSimulation<ProbabilityType>::perform_single_step (void)
{
    const ProbabilityType probability_ratio = walk->compute_probability_ratio_of_random_transition(*rng);

    // phrasing the if statement this way also bails out if
    // probability_ratio == nan
    if (!(probability_ratio >= 0)) {
        // there's generally no good reason to continue once this occurs, but
        // we should at least restore things to a consistent state before
        // throwing an exception
        walk->reject_transition();

        throw invalid_probability_error<ProbabilityType>(probability_ratio);
    }

    ++m_steps;
#if defined(VMC_METROPOLIS_SIMULATION_LOGGING)
    if (m_steps % 200 == 0)
        std::cerr << m_steps << " steps complete" << std::endl;
#endif

    if (probability_ratio >= 1 || probability_ratio > rng->random_uniform01()) {
        // accept transition
#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
        std::cerr << "A" << std::endl;
#endif
        walk->accept_transition();
        ++m_steps_accepted;
        return true;
    } else {
        // reject transition
#if defined(DEBUG_VMC_METROPOLIS_SIMULATION) || defined(DEBUG_VMC_ALL)
        std::cerr << "-" << std::endl;
#endif
        walk->reject_transition();
        if (probability_ratio == 0)
            ++m_steps_fully_rejected;
        return false;
    }
}

template <typename ProbabilityType>
void MetropolisSimulation<ProbabilityType>::perform_initialization (unsigned int initialization_sweeps)
{
    // do initialization sweeps
    for (unsigned int i = 0; i < initialization_sweeps; ++i)
        perform_single_step();
    // initialize the measurements
    for (auto m = measurements.begin(); m != measurements.end(); ++m)
        (*m)->initialize(*walk);
}

#define VMC_SUPPORTED_TYPE(amplitude_type) template class MetropolisSimulation<typename RealPart<amplitude_type>::type>
#include "vmc-supported-types.hpp"
