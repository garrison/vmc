#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>
#include <boost/random.hpp>

typedef boost::mt19937 rng_class;
typedef unsigned long long rng_seed_t;

typedef double probability_t;

typedef double real_amplitude_t;
typedef std::complex<double> complex_amplitude_t;

class accumulator_t
{
public:
    accumulator_t (double x=0)
	: accum(x)
	{
	}

    operator double () const
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
