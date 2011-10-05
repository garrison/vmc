#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>

typedef double probability_t;
typedef std::complex<double> amplitude_t;

class accumulator_t
{
public:
    accumulator_t (double x=0)
	: accum(x)
	{
	}

    operator double ()
	{
	    return accum;
	}

    accumulator_t & operator+= (double x)
	{
	    accum += x;
	    return *this;
	}

    accumulator_t & operator-= (double x)
	{
	    return operator+=(-x);
	}

private:
    long double accum;
};

#endif
