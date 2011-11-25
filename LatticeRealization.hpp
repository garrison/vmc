#ifndef _LATTICE_REALIZATION_HPP
#define _LATTICE_REALIZATION_HPP

#include <cmath>

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/rational.hpp>
#include <boost/math/constants/constants.hpp>

#include "vmc-typedefs.hpp"
#include "NDLattice.hpp"
#include "safe-modulus.hpp"

// this class currently assumes a bravais lattice.  if it's not a bravais
// lattice we will also want to store vectors that point to each site within a
// unit cell.

template<std::size_t DIM>
class LatticeRealization : public NDLattice<DIM>
{
public:
    LatticeRealization (const boost::array<int, DIM> &length_, const boost::array<boost::array<real_position_t, DIM>, DIM> &primitive_vectors_, const boost::array<boost::array<real_position_t, DIM>, DIM> &reciprocal_primitive_vectors_)
        : NDLattice<DIM>(length_, 1),
          primitive_vectors(primitive_vectors_),
          reciprocal_primitive_vectors(reciprocal_primitive_vectors_)
        {
#ifndef BOOST_DISABLE_ASSERTS
            // check that reciprocal and primitive vectors make sense wrt orthogonality condition
            const real_position_t two_pi = 2 * boost::math::constants::pi<real_position_t>();
            for (unsigned int i = 0; i < DIM; ++i) {
                for (unsigned int j = 0; j < DIM; ++j) {
                    real_position_t dot_productish = (i == j) ? -two_pi : 0;
                    for (unsigned int k = 0; k < DIM; ++k)
                        dot_productish += primitive_vectors[i][k] * reciprocal_primitive_vectors[j][k];
                    BOOST_ASSERT(std::abs(dot_productish < .000001));
                }
            }
#endif
        }

    void map_momentum_to_brillouin_zone (boost::array<boost::rational<int>, DIM> &momentum) const
        {
            // first, get everything in interval (-1/2, 1/2]
            for (unsigned int i = 0; i < DIM; ++i) {
                int numerator = momentum[i].numerator(), denominator = momentum[i].denominator();
                do_safe_modulus(numerator, denominator);
                momentum[i] = boost::rational<int>(numerator, denominator);
                if (momentum[i] > boost::rational<int>(1, 2))
                    momentum[i] -= 1;
                BOOST_ASSERT(momentum[i] > boost::rational<int>(-1, 2) && momentum[i] <= boost::rational<int>(1, 2));
            }

            // FIXME: now do a greedy search to get as close to the origin as possible
        }

    const boost::array<boost::array<real_position_t, DIM>, DIM> primitive_vectors, reciprocal_primitive_vectors;
};

#endif
