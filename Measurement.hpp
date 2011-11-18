#ifndef _MEASUREMENT_HPP
#define _MEASUREMENT_HPP

#include <boost/assert.hpp>
#include <boost/noncopyable.hpp>

template <class Walk_T>
class Measurement : boost::noncopyable
// abstract base class
{
public:
    virtual ~Measurement (void)
        {
        }

    void initialize (const Walk_T &walk)
        {
            BOOST_ASSERT(!initialized);
            initialize_(walk);
#ifndef BOOST_DISABLE_ASSERTS
            initialized = true;
#endif
        }

    void measure (const Walk_T &walk)
        {
            BOOST_ASSERT(initialized);
            BOOST_ASSERT(!measurement_in_progress);
            ++m_measurements_completed;
#ifndef BOOST_DISABLE_ASSERTS
            measurement_in_progress = true;
#endif
            measure_(walk);
#ifndef BOOST_DISABLE_ASSERTS
            measured = true;
            measurement_in_progress = false;
#endif
        }

    void repeat_measurement (const Walk_T &walk)
        {
            BOOST_ASSERT(initialized && measured);
            BOOST_ASSERT(!measurement_in_progress);
            ++m_measurements_completed;
#ifndef BOOST_DISABLE_ASSERTS
            measurement_in_progress = true;
#endif
            repeat_measurement_(walk);
#ifndef BOOST_DISABLE_ASSERTS
            measurement_in_progress = false;
#endif
        }

protected:
    Measurement (void)
        : m_measurements_completed(0)
#ifndef BOOST_DISABLE_ASSERTS
        , initialized(false)
        , measured(false)
        , measurement_in_progress(false)
#endif
        {
        }

    unsigned int get_measurements_completed (void) const
        {
            return m_measurements_completed;
        }

private:
    virtual void initialize_ (const Walk_T &walk)
        {
            // by default, do nothing.
            (void) walk;
        }

    virtual void measure_ (const Walk_T &walk) = 0;

    virtual void repeat_measurement_ (const Walk_T &walk)
        {
            // by default, call measure_(), NOT measure().  Calling measure()
            // would increment m_measurements_completed redundantly.
            measure_(walk);
        }

    unsigned int m_measurements_completed;

#ifndef BOOST_DISABLE_ASSERTS
    bool initialized, measured, measurement_in_progress;
#endif
};

#endif
