#ifndef _LATTICE_REALIZATION_HPP
#define _LATTICE_REALIZATION_HPP

#include <cmath>

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/rational.hpp>
#include <boost/math/constants/constants.hpp>

#include "vmc-typedefs.hpp"
#include "Lattice.hpp"
#include "vmc-math-utils.hpp"

// this class currently assumes a bravais lattice.  if it's not a bravais
// lattice we will also want to store vectors that point to each site within a
// unit cell.

/**
 * A lattice that exists with definite primitive vectors in real space.
 *
 * Most calculations only depend on what dimension the lattice is, but
 * occasionally we want to represent a lattice with sites in real space.  For
 * instance, only if we know the primitive vectors can we determine the
 * magnitude of each momenta, so having a LatticeRealization is necessary if we
 * e.g. want to fill the M lowest momenta.
 */
template<std::size_t DIM>
class LatticeRealization : public NDLattice<DIM>
{
public:
    /**
     * Constructor.
     *
     * @param length_ size of lattice
     *
     * @param primitive_vectors_ primitive vectors of lattice
     *
     * @param reciprocal_primitive_vectors_ primitive vectors of reciprocal
     * lattice, must satisfy \f$ \vec{a}_i \vec{b}_j = 2\pi\delta_{ij}\f$ where
     * \f$\vec{a}_i\f$ are the direct lattice primitive vectors and
     * \f$\vec{b}_j\f$ are the reciprocal lattice primitive vectors
     */
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

    /**
     * Modifies the given momentum site so that it points to the same momentum
     * value but within the first Brillouin zone.
     *
     * FIXME: currently this works only for hypercubic lattices
     */
    void map_momentum_to_brillouin_zone (boost::array<boost::rational<int>, DIM> &momentum) const
        {
            // first, get everything in interval [-1/2, 1/2)
            for (unsigned int i = 0; i < DIM; ++i) {
                int numerator = momentum[i].numerator(), denominator = momentum[i].denominator();
                do_safe_modulus(numerator, denominator);
                momentum[i] = boost::rational<int>(numerator, denominator);
                if (momentum[i] >= boost::rational<int>(1, 2))
                    momentum[i] -= 1;
                BOOST_ASSERT(momentum[i] >= boost::rational<int>(-1, 2) && momentum[i] < boost::rational<int>(1, 2));
            }

            // FIXME: now do a greedy search to get as close to the origin as possible
        }

    const boost::array<boost::array<real_position_t, DIM>, DIM> primitive_vectors, reciprocal_primitive_vectors;
};

#endif
