#ifndef _MEASUREMENT_HPP
#define _MEASUREMENT_HPP

#include <boost/cast.hpp>
#include <boost/assert.hpp>
#include <boost/noncopyable.hpp>

#include "Walk.hpp"

class BaseMeasurement : boost::noncopyable
// common (non-templated) abstract base class
{
public:
    virtual ~BaseMeasurement (void)
        {
        }

    virtual void initialize (const Walk &walk) = 0;

    virtual void step_advanced (const Walk &walk) = 0;

    virtual void step_repeated (const Walk &walk) = 0;

    virtual bool is_valid_walk (const Walk &walk) = 0;
};

template <class Walk_T>
class Measurement : public BaseMeasurement
// abstract base class
{
public:
    /**
     * Initializes the object for taking measurements.
     *
     * This is done once the walk is brought to "equilibrium."  After this is
     * called, either measure() or repeat_measurement() should be called after
     * every step.
     */
    void initialize (const Walk &walk)
        {
            BOOST_ASSERT(!initialized);
            initialize_(*boost::polymorphic_downcast<const Walk_T *>(&walk));
            initialized = true;
        }

    /**
     * Once we wish to begin taking measurements, this should be called any
     * time a step is taken that results in an actual move in Monte Carlo
     * space.
     */
    void step_advanced (const Walk &walk)
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
    void step_repeated (const Walk &walk)
        {
            BOOST_ASSERT(first_step_has_been_completed);
            step_completed(walk);
        }

    bool is_valid_walk (const Walk &walk)
        {
            return bool(dynamic_cast<const Walk_T *>(&walk));
        }

protected:
    Measurement (unsigned int steps_per_measurement)
        : m_steps_per_measurement(steps_per_measurement),
          m_steps_since_last_measurement(0),
          m_state_changed_since_last_measurement(false),
          initialized(false),
          first_step_has_been_completed(false),
          measurement_in_progress(false)
        {
            BOOST_ASSERT(steps_per_measurement > 0);
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

    /**
     * Returns true if the measurement has been initialized with a Walk
     *
     * @see initialize()
     */
    bool is_initialized (void) const
        {
            return initialized;
        }

private:
    /**
     * This is called any time a step is made.  An actual measurement will be
     * performed here if one is due.
     */
    void step_completed (const Walk &walk_)
        {
            BOOST_ASSERT(initialized);
            BOOST_ASSERT(!measurement_in_progress);
            const Walk_T &walk = *boost::polymorphic_downcast<const Walk_T *>(&walk_);
            ++m_steps_since_last_measurement;
            if (m_steps_since_last_measurement % m_steps_per_measurement == 0) {
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
    unsigned int m_steps_since_last_measurement;
    bool m_state_changed_since_last_measurement;

    bool initialized;

    // these are used only for assertions
    bool first_step_has_been_completed, measurement_in_progress;
};

#endif
