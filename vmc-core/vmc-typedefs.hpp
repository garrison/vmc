#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>

typedef double real_t;
typedef std::complex<real_t> complex_t;

typedef real_t probability_t;

typedef complex_t amplitude_t;

static const unsigned int MAX_DIMENSION = 3;
static const unsigned int MAX_MOVE_SIZE = 4;

// the below code is from <http://stackoverflow.com/a/17902439/1558890>.  Since
// std::make_unique exists only in c++1y (and later), we include this here for
// c++11 compatibility.

#if __cplusplus == 201103L

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace std {
    template<class T> struct _Unique_if {
        typedef unique_ptr<T> _Single_object;
    };

    template<class T> struct _Unique_if<T[]> {
        typedef unique_ptr<T[]> _Unknown_bound;
    };

    template<class T, size_t N> struct _Unique_if<T[N]> {
        typedef void _Known_bound;
    };

    template<class T, class... Args>
        typename _Unique_if<T>::_Single_object
        make_unique(Args&&... args) {
            return unique_ptr<T>(new T(std::forward<Args>(args)...));
        }

    template<class T>
        typename _Unique_if<T>::_Unknown_bound
        make_unique(size_t n) {
            typedef typename remove_extent<T>::type U;
            return unique_ptr<T>(new U[n]());
        }

    template<class T, class... Args>
        typename _Unique_if<T>::_Known_bound
        make_unique(Args&&...) = delete;
}

#endif

#endif
