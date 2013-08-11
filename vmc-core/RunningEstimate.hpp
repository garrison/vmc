#ifndef _VMC_RUNNING_ESTIMATE_HPP
#define _VMC_RUNNING_ESTIMATE_HPP

#ifndef BOOST_NUMERIC_FUNCTIONAL_STD_COMPLEX_SUPPORT
#define BOOST_NUMERIC_FUNCTIONAL_STD_COMPLEX_SUPPORT
#endif

#include <cassert>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/count.hpp>

#include "vmc-typedefs.hpp"

template <typename T>
class RunningEstimate
{
public:
    typedef typename boost::numeric::functional::average<T, std::size_t>::result_type result_t;

    virtual ~RunningEstimate (void)
        {
        }

    virtual void add_value (T value)
        {
            m_cumulative_acc(value);
        }

    /**
     * Returns the average of all measurements
     */
    result_t get_cumulative_result (void) const
        {
            assert(boost::accumulators::count(m_cumulative_acc) > 0);
            return boost::accumulators::mean(m_cumulative_acc);
        }

    /**
     * Returns the cumulative number of samples in the history of this
     * estimator
     */
    unsigned int get_num_cumulative_values (void) const
        {
            return boost::accumulators::count(m_cumulative_acc);
        }

protected:
    T get_cumulative_total_value (void) const
        {
            return boost::accumulators::sum(m_cumulative_acc);
        }

private:
    typedef boost::accumulators::accumulator_set<T, boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::sum, boost::accumulators::tag::count> > accumulator_t;

    accumulator_t m_cumulative_acc;
};

#endif
