#ifndef _ARRAY_UTIL_HPP
#define _ARRAY_UTIL_HPP

#include <cstddef>

#include <boost/array.hpp>

// utility functions for easily creating a small array

template<typename T, std::size_t N>
inline boost::array<T, N> make_array (const T &v)
{
    boost::array<T, N> rv;
    rv.assign(v);
    return rv;
}

template<typename T>
inline boost::array<T, 1> make_array (const T &v1)
{
    const boost::array<T, 1> rv = { { v1 } };
    return rv;
}

template<typename T>
inline boost::array<T, 2> make_array (const T &v1, const T &v2)
{
    const boost::array<T, 2> rv = { { v1, v2 } };
    return rv;
}

template<typename T>
inline boost::array<T, 3> make_array (const T &v1, const T &v2, const T &v3)
{
    const boost::array<T, 3> rv = { { v1, v2, v3 } };
    return rv;
}

#endif
