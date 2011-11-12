#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>
#include <boost/random.hpp>

typedef boost::mt19937 rng_class;
typedef unsigned long long rng_seed_t;

typedef long double real_t;
typedef std::complex<real_t> complex_t;

typedef real_t probability_t;

typedef complex_t amplitude_t;

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

#if defined(__GNUC__)

#include <Eigen/Core>

namespace Eigen {

template<> struct NumTraits<__float128>
{
    typedef __float128 Real;
    typedef __float128 NonInteger;
    typedef __float128 Nested;

    enum {
	IsComplex = 0,
	IsInteger = 0,
	IsSigned,
	ReadCost = 1,
	AddCost = 1,
	MulCost = 1,
	RequireInitialization = 0
    };
};

}

// FIXME: overloading std namespace is evil ...
namespace std {

inline __float128 abs (const __float128 &v)
{
    return std::abs(static_cast<double>(v));
}

inline __float128 abs (const std::complex<__float128> &v)
{
    return std::abs(std::complex<double>(v.real(), v.imag()));
}

inline __float128 arg (const std::complex<__float128> &v)
{
    return std::arg(std::complex<double>(v.real(), v.imag()));
}

inline __float128 atan2 (const __float128 &v1, const __float128 &v2)
{
    return std::atan2(static_cast<double>(v1), static_cast<double>(v2));
}

}

#endif

#endif
