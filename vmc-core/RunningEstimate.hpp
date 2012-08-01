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
    // get_result()
    template <bool, class>
    class ResultType;

    // if we are estimating an integral type, get_result() should return a
    // real_t
    template <class T1>
    class ResultType<true, T1>
    {
    public:
        typedef real_t type;
    };

    // if we are estimating a non-integral type, get_result() should return
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
        : m_total_value(T(0)),
          m_num_values(0)
        {
        }

    virtual ~RunningEstimate (void)
        {
        }

    void add_value (T value)
        {
            m_total_value += value;
            ++m_num_values;
        }

    result_t get_result (void) const
        {
            BOOST_ASSERT(m_num_values > 0);
            return m_total_value / real_t(m_num_values);
        }

    unsigned int get_num_values (void) const
        {
            return m_num_values;
        }

protected:
    // fixme: will we actually use this function?
    T get_total_value (void) const
        {
            return m_total_value;
        }

private:
    T m_total_value;
    unsigned int m_num_values;
};

#endif
