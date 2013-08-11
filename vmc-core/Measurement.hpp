#ifndef _VMC_MEASUREMENT_HPP
#define _VMC_MEASUREMENT_HPP

#include <cassert>

#include <boost/cast.hpp>
#include <boost/noncopyable.hpp>

#include "Walk.hpp"

template <typename ProbabilityType>
class BaseMeasurement : boost::noncopyable
// common (non-templated) abstract base class
{
public:
    typedef Walk<ProbabilityType> BaseWalkType;

    virtual ~BaseMeasurement (void)
        {
        }

    virtual void initialize (const BaseWalkType &walk) = 0;

    virtual void step_advanced (const BaseWalkType &walk) = 0;

    virtual void step_repeated (const BaseWalkType &walk) = 0;

    virtual bool is_valid_walk (const BaseWalkType &walk) = 0;
};

template <class _WalkType>
class Measurement : public BaseMeasurement<typename _WalkType::ProbabilityType>
// abstract base class
{
public:
    typedef _WalkType WalkType;
    typedef typename WalkType::ProbabilityType ProbabilityType;
    typedef Walk<ProbabilityType> BaseWalkType;

    /**
     * Initializes the object for taking measurements.
     *
     * This is done once the walk is brought to "equilibrium."  After this is
     * called, either measure() or repeat_measurement() should be called after
     * every step.
     */
    virtual void initialize (const BaseWalkType &walk) override final
        {
            assert(!initialized);
            assert(this->is_valid_walk(walk));
            initialize_(*boost::polymorphic_downcast<const WalkType *>(&walk));
            initialized = true;
        }

    /**
     * Once we wish to begin taking measurements, this should be called any
     * time a step is taken that results in an actual move in Monte Carlo
     * space.
     */
    virtual void step_advanced (const BaseWalkType &walk) override final
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
    virtual void step_repeated (const BaseWalkType &walk) override final
        {
            assert(first_step_has_been_completed);
            step_completed(walk);
        }

    virtual bool is_valid_walk (const BaseWalkType &walk) override final
        {
            const WalkType *walkptr = dynamic_cast<const WalkType *>(&walk);
            return walkptr && is_valid_walk_(*walkptr);
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
            assert(steps_per_measurement > 0);
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
    void step_completed (const BaseWalkType &walk_)
        {
            assert(initialized);
            assert(!measurement_in_progress);
            const WalkType &walk = *boost::polymorphic_downcast<const WalkType *>(&walk_);
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

    virtual void initialize_ (const WalkType &walk)
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
    virtual void measure_ (const WalkType &walk) = 0;

    /**
     * This gets called any time we actually want to tally a measurement
     * (i.e. we have already advanced by steps_per_measurement steps).
     *
     * Specifically, this is called when the walk has not actually advanced
     * since the most recent measurement.
     *
     * @see measure_()
     */
    virtual void repeat_measurement_ (const WalkType &walk)
        {
            // by default, simply call measure_()
            measure_(walk);
        }

    virtual bool is_valid_walk_ (const WalkType &walk)
        {
            // no additional constraints by default
            (void) walk;
            return true;
        }

    const unsigned int m_steps_per_measurement;
    unsigned int m_steps_since_last_measurement;
    bool m_state_changed_since_last_measurement;

    bool initialized;

    // these are used only for assertions
    bool first_step_has_been_completed, measurement_in_progress;
};

#endif
