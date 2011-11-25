#ifndef _SAFE_MODULUS_HPP
#define _SAFE_MODULUS_HPP

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

#endif
