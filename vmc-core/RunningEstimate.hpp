#ifndef _RUNNING_ESTIMATE_HPP
#define _RUNNING_ESTIMATE_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/assert.hpp>

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

    RunningEstimate (void)
        : m_recent_total_value(T(0)),
          m_cumulative_total_value(T(0)),
          m_num_recent_values(0),
          m_num_cumulative_values(0)
        {
        }

    virtual ~RunningEstimate (void)
        {
        }

    virtual void add_value (T value)
        {
            m_recent_total_value += value;
            ++m_num_recent_values;

            m_cumulative_total_value += value;
            ++m_num_cumulative_values;
        }

    /**
     * Returns the average of all measurements since the most recent reset
     */
    result_t get_recent_result (void) const
        {
            BOOST_ASSERT(m_num_recent_values > 0);
            return m_recent_total_value / real_t(m_num_recent_values);
        }

    /**
     * Returns the average of all measurements, regardless of whether the
     * simulation has been reset
     */
    result_t get_cumulative_result (void) const
        {
            BOOST_ASSERT(m_num_cumulative_values > 0);
            return m_cumulative_total_value / real_t(m_num_cumulative_values);
        }

    /**
     * Returns the number of samples since the last reset
     */
    unsigned int get_num_recent_values (void) const
        {
            return m_num_recent_values;
        }

    /**
     * Returns the cumulative number of samples in the history of this
     * estimator
     */
    unsigned int get_num_cumulative_values (void) const
        {
            return m_num_cumulative_values;
        }

    /**
     * Resets the estimator
     */
    void reset (void)
        {
            m_num_recent_values = 0;
            m_recent_total_value = T(0);
        }

protected:
    T get_cumulative_total_value (void) const
        {
            return m_cumulative_total_value;
        }

private:
    T m_recent_total_value, m_cumulative_total_value;
    unsigned int m_num_recent_values, m_num_cumulative_values;
};

#endif
