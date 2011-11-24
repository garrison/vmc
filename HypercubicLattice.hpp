#ifndef _HYPERCUBIC_LATTICE_HPP
#define _HYPERCUBIC_LATTICE_HPP

#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>

#include "LatticeRealization.hpp"
#include "vmc-typedefs.hpp"

template<std::size_t DIM>
class HypercubicLattice : public LatticeRealization<DIM>
{
public:
    HypercubicLattice (const boost::array<int, DIM> &length_, real_position_t a=1)
        : LatticeRealization<DIM>(length_, generate_vectors(a), generate_vectors(2 * boost::math::constants::pi<real_position_t>() / a))
        {
        }

private:
    static boost::array<boost::array<real_position_t, DIM>, DIM> generate_vectors(real_position_t vlength)
        {
            boost::array<boost::array<real_position_t, DIM>, DIM> rv;
            for (unsigned int i = 0; i < DIM; ++i) {
                rv[i].assign(0);
                rv[i][i] = vlength;
            }
            return rv;
        }
};

#endif
