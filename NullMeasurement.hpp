#ifndef _NULL_MEASUREMENT_HPP
#define _NULL_MEASUREMENT_HPP

#include "Measurement.hpp"

template <class Walk_T>
class NullMeasurement : public Measurement<Walk_T>
{
public:
    typedef int measurement_value_t;

    NullMeasurement (void)
	{
	}

    measurement_value_t get (void) const
	{
	    return 0;
	}

private:
    void measure_ (const Walk_T &walk)
	{
	    (void) walk;
	}
};

#endif
