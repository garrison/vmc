#ifndef _VMC_TYPEDEFS_HPP
#define _VMC_TYPEDEFS_HPP

#include <complex>

typedef double real_t;
typedef std::complex<real_t> complex_t;

typedef real_t probability_t;

typedef complex_t amplitude_t;
typedef amplitude_t phase_t;

static const unsigned int MAX_DIMENSION = 3;
static const unsigned int MAX_MOVE_SIZE = 4;

#endif
