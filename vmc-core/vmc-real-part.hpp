#ifndef _VMC_REAL_PART_HPP
#define _VMC_REAL_PART_HPP

#include <complex>

// see http://stackoverflow.com/q/14925544/1558890

template <typename T>
struct RealPart
{
    typedef T type;
};

template <typename T>
struct RealPart<std::complex<T> >
{
    typedef T type;
};

#endif
