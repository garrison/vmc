#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>

typedef double real_t;
typedef std::complex<real_t> complex_t;

typedef real_t probability_t;
typedef real_t real_position_t;

typedef complex_t amplitude_t;
typedef complex_t phase_t;

static const unsigned int MAX_DIMENSION = 3;
static const unsigned int MAX_MOVE_SIZE = 2;

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
