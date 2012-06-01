#include <boost/assert.hpp>

#include "allowed-momentum.hpp"
#include "Lattice.hpp"
#include "vmc-math-utils.hpp"

lw_vector<boost::rational<int>, MAX_DIMENSION> allowed_momentum (const lw_vector<int, MAX_DIMENSION> &momentum_site, const Lattice &lattice, const BoundaryConditions &bcs)
{
    BOOST_ASSERT(lattice.n_dimensions() == bcs.size());
    BOOST_ASSERT(lattice.n_dimensions() == momentum_site.size());
    lw_vector<boost::rational<int>, MAX_DIMENSION> rv(lattice.n_dimensions());
    for (unsigned int i = 0; i < lattice.n_dimensions(); ++i) {
        BOOST_ASSERT(momentum_site[i] >= 0 && momentum_site[i] < lattice.dimensions[i]);
        rv[i] = (bcs[i].p_floor() + momentum_site[i]) * boost::rational<int>(1, lattice.dimensions[i]);
    }
    return rv;
}

#if 0
lw_vector<boost::rational<int>, MAX_DIMENSION> allowed_momentum_safe (const lw_vector<int, MAX_DIMENSION> &momentum_site, const Lattice &lattice, const BoundaryConditions &bcs)
{
    BOOST_ASSERT(lattice.n_dimensions() == bcs.size());
    BOOST_ASSERT(lattice.n_dimensions() == momentum_site.size());
    lw_vector<int, MAX_DIMENSION> momentum_site_(momentum_site);
    for (unsigned int i = 0; i < lattice.n_dimensions(); ++i)
        do_safe_modulus(momentum_site[i], lattice.dimensions[i]);
    return allowed_momentum(momentum_site_, lattice, bcs);
}
#endif
