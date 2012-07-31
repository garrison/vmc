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
     * Once we wish to begin taking measurements, this should be called any
     * time a step is taken that results in an actual move in Monte Carlo
     * space.
     */
    void step_advanced (const Walk_T &walk)
        {
            first_step_has_been_completed = true;
            m_state_changed_since_last_measurement = true;
            step_completed(walk);
        }

    /**
     * Once we wish to begin taking measurements, this should be called any
     * time a step is taken in which the current state of the walk is identical
     * to the previous state.
     */
    void step_repeated (const Walk_T &walk)
        {
            BOOST_ASSERT(first_step_has_been_completed);
            step_completed(walk);
        }

protected:
    Measurement (unsigned int steps_per_measurement)
        : m_steps_per_measurement(steps_per_measurement),
          m_measurements_completed(0),
          m_steps_since_last_measurement(0),
          m_state_changed_since_last_measurement(false),
          initialized(false),
          first_step_has_been_completed(false),
          measurement_in_progress(false)
        {
            BOOST_ASSERT(steps_per_measurement > 0);
        }

    /**
     * Returns how many measurements have been tallied
     */
    unsigned int get_measurements_completed (void) const
        {
            return m_measurements_completed;
        }

    /**
     * Returns the steps that will be taken between each measurement.
     *
     * This is constant once a Measurement is instantiated.
     */
    unsigned int get_steps_per_measurement (void) const
        {
            return m_steps_per_measurement;
        }

    /**
     * Returns the number of steps since the last measurement
     */
    unsigned int get_steps_since_last_measurement (void) const
        {
            return m_steps_since_last_measurement;
        }

private:
    /**
     * This is called any time a step is made.  An actual measurement will be
     * performed here if one is due.
     */
    void step_completed (const Walk_T &walk)
        {
            BOOST_ASSERT(initialized);
            BOOST_ASSERT(!measurement_in_progress);
            ++m_steps_since_last_measurement;
            if (m_steps_since_last_measurement % m_steps_per_measurement == 0) {
                ++m_measurements_completed;
                measurement_in_progress = true;
                if (m_state_changed_since_last_measurement) {
                    measure_(walk);
                } else {
                    repeat_measurement_(walk);
                }
                measurement_in_progress = false;
                m_steps_since_last_measurement = 0;
                m_state_changed_since_last_measurement = false;
            }
        }

    virtual void initialize_ (const Walk_T &walk)
        {
            // by default, do nothing.
            (void) walk;
        }

    /**
     * This gets called any time we actually want to tally a measurement
     * (i.e. we have already advanced by steps_per_measurement steps).
     *
     * Specifically, this is called when the walk has actually advanced since
     * the most recent measurement.
     *
     * @see repeat_measurement_()
     */
    virtual void measure_ (const Walk_T &walk) = 0;

    /**
     * This gets called any time we actually want to tally a measurement
     * (i.e. we have already advanced by steps_per_measurement steps).
     *
     * Specifically, this is called when the walk has not actually advanced
     * since the most recent measurement.
     *
     * @see measure_()
     */
    virtual void repeat_measurement_ (const Walk_T &walk)
        {
            // by default, simply call measure_()
            measure_(walk);
        }

    const unsigned int m_steps_per_measurement;
    unsigned int m_measurements_completed, m_steps_since_last_measurement;
    bool m_state_changed_since_last_measurement;

    // these are used only for assertions
    bool initialized, first_step_has_been_completed, measurement_in_progress;
};

#endif
