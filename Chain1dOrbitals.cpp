#include <cmath>

#include <boost/assert.hpp>
#include <boost/math/constants/constants.hpp>

#include "Chain1dOrbitals.hpp"
#include "HypercubicLattice.hpp"

bool Chain1dOrbitals::lattice_makes_sense (const Lattice &lattice) const
{
    return bool(dynamic_cast<const HypercubicLattice<1> *>(&lattice));
}

amplitude_t Chain1dOrbitals::calculate_phi (unsigned int n, unsigned int r, const Lattice &lattice) const
{
    BOOST_ASSERT(this->lattice_makes_sense(lattice));
    const real_t two_pi = 2 * boost::math::constants::pi<real_t>();
    int kbar = ((n + 1) / 2) * ((n & 1) * -2 + 1); // fill each k value in order, alternating +k, -k
    const complex_t im_unit(0, 1);
    return std::exp(im_unit * complex_t(two_pi * kbar * r / lattice.total_sites()));
}
