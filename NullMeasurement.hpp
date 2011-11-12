#ifndef NULL_MEASUREMENT_HPP
#define NULL_MEASUREMENT_HPP

template <class Walk_T>
class NullMeasurement
{
public:
    typedef int measurement_value_t;

    NullMeasurement (const Walk_T &walk)
	{
	    (void) walk;
	}

    void measure (const Walk_T &walk)
	{
	    (void) walk;
	}

    measurement_value_t get (unsigned int measurements_completed) const
	{
	    (void) measurements_completed;
	    return 0;
	}
};

#endif
