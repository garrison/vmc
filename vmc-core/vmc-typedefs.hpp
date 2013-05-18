#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>

#ifdef VMC_REAL_T
typedef VMC_REAL_T real_t;
#else
typedef double real_t;
#endif

typedef std::complex<real_t> complex_t;

typedef real_t probability_t;

#ifdef VMC_AMPLITUDE_T
typedef VMC_AMPLITUDE_T amplitude_t;
#else
typedef complex_t amplitude_t;
#endif

typedef amplitude_t phase_t;

static const unsigned int MAX_DIMENSION = 3;
static const unsigned int MAX_MOVE_SIZE = 4;

#endif
