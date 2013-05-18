#ifndef _VMC_MATH_UTILS_HPP
#define _VMC_MATH_UTILS_HPP

#include <cmath>

#include <boost/assert.hpp>

#include "vmc-typedefs.hpp"
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

/**
 * Same as std::pow(), but returns zero instead of domain error for 0 ^ z when z < 0
 */
static inline real_t safe_pow (const real_t &base, const real_t &exponent)
{
    if (base == 0.0)
        return 0.0;
    return std::pow(base, exponent);
}

/**
 * Returns base * abs(base) ^ (exponent - 1)
 */
template <typename T>
static inline T complex_pow (const T &base, const real_t &exponent)
{
    T rv(base);
    if (exponent != 1.0) {
        rv *= safe_pow(std::abs(base), exponent - 1.0);
    }
    return rv;
}

template <typename T>
static inline Big<T> complex_pow (const Big<T> &base, const real_t &exponent)
{
    return Big<T>(complex_pow(base.get_base(), exponent), base.get_exponent() * exponent);
}

#endif
