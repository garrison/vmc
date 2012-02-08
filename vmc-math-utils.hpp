#ifndef _VMC_MATH_UTILS_HPP
#define _VMC_MATH_UTILS_HPP

#include <cmath>

#include "vmc-typedefs.hpp"

// the modulus of a negative number is not consistently defined in C/C++, so
// this function is defined to give the behavior one would logically expect

static inline void do_safe_modulus (int &n, int d)
{
    while (n < d)
        n += d;
    n %= d;
}

static inline int safe_modulus (int n, int d)
{
    do_safe_modulus(n, d);
    return n;
}

/**
 * Returns base * abs(base) ^ (exponent - 1)
 */
static inline complex_t complex_pow (complex_t base, real_t exponent)
{
    complex_t rv = base;
    if (exponent != 1.0)
        rv *= std::pow(std::abs(base), exponent - 1.0);
    return rv;
}

#endif
