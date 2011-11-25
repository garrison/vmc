#ifndef _ALLOWED_MOMENTUM_HPP
#define _ALLOWED_MOMENTUM_HPP

#include <cstddef>

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/rational.hpp>

#include "NDLattice.hpp"
#include "safe-modulus.hpp"

// returns things in terms of the primitive vectors of the reciprocal lattice
template<std::size_t DIM>
boost::array<boost::rational<int>, DIM> allowed_momentum (const boost::array<int, DIM> &momentum_site, const NDLattice<DIM> &lattice, const typename NDLattice<DIM>::BoundaryConditions &bcs)
{
    boost::array<boost::rational<int>, DIM> rv;
    for (unsigned int i = 0; i < DIM; ++i) {
        BOOST_ASSERT(momentum_site[i] >= 0 && momentum_site[i] < lattice.length[i]);
        rv[i] = (bcs[i].p() + momentum_site[i] - 1) * boost::rational<int>(1, lattice.length[i]);
    }
    return rv;
}

#if 0
// same thing, but performs any necessary wrapping in the integer indices
template<std::size_t DIM>
boost::array<boost::rational<int>, DIM> allowed_momentum_safe (const boost::array<int, DIM> &momentum_site, const NDLattice<DIM> &lattice, const typename NDLattice<DIM>::BoundaryConditions &bcs)
{
    boost::array<int, DIM> momentum_site_(momentum_site);
    for (unsigned int i = 0; i < DIM; ++i)
        do_safe_modulus(momentum_site[i], lattice.length[i]);
    return allowed_momentum(momentum_site_, lattice, bcs);
}
#endif

#endif
