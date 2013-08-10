#ifndef _VMC_MATH_UTILS_HPP
#define _VMC_MATH_UTILS_HPP

#include <cmath>
#include <complex>

#include <boost/assert.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_complex.hpp>

#include "vmc-real-part.hpp"
#include "Big.hpp"

// the modulus of a negative number is not consistently defined in C/C++, so
// this function is defined to give the behavior one would logically expect

// from http://stackoverflow.com/a/4003287/1558890
static inline void do_safe_modulus (int &n, int d)
{
    BOOST_ASSERT(d > 0);
    n %= d;
    if (n < 0)
        n += d;
}

static inline int safe_modulus (int n, int d)
{
    do_safe_modulus(n, d);
    return n;
}

template <typename T>
static inline T safe_fmod (const T &n, const T &d)
{
    BOOST_ASSERT(d > 0);
    using std::fmod;
    T rv = fmod(n, d);
    if (rv < 0)
        rv += d;
    return rv;
}

/**
 * Same as std::pow(), but returns zero instead of domain error for 0 ^ z when z < 0
 */
template <typename RealType>
static inline typename boost::disable_if<boost::is_complex<RealType>, RealType>::type
safe_pow (const RealType &base, const RealType &exponent)
{
    if (base == 0.0)
        return 0.0;
    using std::pow;
    return pow(base, exponent);
}

/**
 * Returns base * abs(base) ^ (exponent - 1)
 */
template <typename T>
static inline T complex_pow (const T &base, const typename RealPart<T>::type &exponent)
{
    using std::abs;
    T rv(base);
    if (exponent != 1.0) {
        rv *= safe_pow(abs(base), exponent - 1.0);
    }
    return rv;
}

template <typename T, typename T2>
static inline typename boost::enable_if<boost::is_convertible<T2, typename RealPart<T>::type>, Big<T> >::type
complex_pow (const Big<T> &base, const T2 &exponent)
{
    return Big<T>(complex_pow(base.get_base(), exponent), base.get_exponent() * exponent);
}

// std::conj() always returns a complex type, but we want a function that
// always returns the type it is given (e.g. returns real if given a real).
// vmc_conj() is meant to fill in this gap.

template <typename T>
static inline T
vmc_conj (const T &v)
{
    return v;
}

template <typename T>
static inline std::complex<T>
vmc_conj (const std::complex<T> &v)
{
    using std::conj;
    return conj(v);
}

#endif
