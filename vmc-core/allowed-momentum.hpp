#ifndef _ALLOWED_MOMENTUM_HPP
#define _ALLOWED_MOMENTUM_HPP

#include <boost/rational.hpp>

#include "lw_vector.hpp"
#include "vmc-typedefs.hpp"
#include "BoundaryCondition.hpp"

class Lattice;

// returns things in terms of the primitive vectors of the reciprocal lattice
extern lw_vector<boost::rational<int>, MAX_DIMENSION> allowed_momentum (const lw_vector<int, MAX_DIMENSION> &momentum_site, const Lattice &lattice, const BoundaryConditions &bcs);

#if 0
// same thing, but performs any necessary wrapping in the integer indices
extern lw_vector<boost::rational<int>, MAX_DIMENSION> allowed_momentum_safe (const lw_vector<int, MAX_DIMENSION> &momentum_site, const Lattice &lattice, const BoundaryConditions &bcs);
#endif

#endif
