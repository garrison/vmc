#ifndef _RUNNING_ESTIMATE_HPP
#define _RUNNING_ESTIMATE_HPP

#ifndef BOOST_NUMERIC_FUNCTIONAL_STD_COMPLEX_SUPPORT
#define BOOST_NUMERIC_FUNCTIONAL_STD_COMPLEX_SUPPORT
#endif

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/assert.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/count.hpp>

#include "vmc-typedefs.hpp"

template <typename T>
class RunningEstimate
{
private:
    // The point of the following section of code is to use C++ "partial
    // template specialization" to determine the correct return type for
    // get_*_result()
    template <bool, class>
    class ResultType;

    // if we are estimating an integral type, get_*_result() should return a
    // real_t
    template <class T1>
    class ResultType<true, T1>
    {
    public:
        typedef real_t type;
    };

    // if we are estimating a non-integral type, get_*_result() should return
    // that type
    template <class T1>
    class ResultType<false, T1>
    {
    public:
        typedef T1 type;
    };

public:
    // define result_t, making use of the template specialization above
    typedef typename ResultType<boost::is_integral<T>::value, T>::type result_t;

    virtual ~RunningEstimate (void)
        {
        }

    virtual void add_value (T value)
        {
            m_recent_acc(value);
            m_cumulative_acc(value);
        }

    /**
     * Returns the average of all measurements since the most recent reset
     */
    result_t get_recent_result (void) const
        {
            BOOST_ASSERT(boost::accumulators::count(m_recent_acc) > 0);
            return boost::accumulators::mean(m_recent_acc);
        }

    /**
     * Returns the average of all measurements, regardless of whether the
     * simulation has been reset
     */
    result_t get_cumulative_result (void) const
        {
            BOOST_ASSERT(boost::accumulators::count(m_cumulative_acc) > 0);
            return boost::accumulators::mean(m_cumulative_acc);
        }

    /**
     * Returns the number of samples since the last reset
     */
    unsigned int get_num_recent_values (void) const
        {
            return boost::accumulators::count(m_recent_acc);
        }

    /**
     * Returns the cumulative number of samples in the history of this
     * estimator
     */
    unsigned int get_num_cumulative_values (void) const
        {
            return boost::accumulators::count(m_cumulative_acc);
        }

    /**
     * Resets the estimator
     */
    void reset (void)
        {
            m_recent_acc = accumulator_t();
        }

protected:
    T get_cumulative_total_value (void) const
        {
            return boost::accumulators::sum(m_cumulative_acc);
        }

private:
    typedef boost::accumulators::accumulator_set<T, boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::sum, boost::accumulators::tag::count> > accumulator_t;

    accumulator_t m_recent_acc, m_cumulative_acc;
};

#endif
