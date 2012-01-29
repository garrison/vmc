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

    /**
     * Initializes the object for taking measurements.
     *
     * This is done once the walk is brought to "equilibrium."  After this is
     * called, either measure() or repeat_measurement() should be called after
     * every step.
     */
    void initialize (const Walk_T &walk)
        {
            BOOST_ASSERT(!initialized);
            initialize_(walk);
            initialized = true;
        }

    /**
     * Calculates and tallies a measurement
     */
    void measure (const Walk_T &walk)
        {
            BOOST_ASSERT(initialized);
            BOOST_ASSERT(!measurement_in_progress);
            ++m_measurements_completed;
            measurement_in_progress = true;
            measure_(walk);
            measured = true;
            measurement_in_progress = false;
        }

    /**
     * Tallies again the most recent measurement
     */
    void repeat_measurement (const Walk_T &walk)
        {
            BOOST_ASSERT(initialized && measured);
            BOOST_ASSERT(!measurement_in_progress);
            ++m_measurements_completed;
            measurement_in_progress = true;
            repeat_measurement_(walk);
            measurement_in_progress = false;
        }

protected:
    Measurement (void)
        : m_measurements_completed(0),
          initialized(false),
          measured(false),
          measurement_in_progress(false)
        {
        }

    /**
     * Returns how many measurements have been tallied
     */
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

    bool initialized, measured, measurement_in_progress;
};

#endif
