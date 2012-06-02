#ifndef _LOWEST_MOMENTA_HPP
#define _LOWEST_MOMENTA_HPP

#include <cstddef>
#include <vector>
#include <algorithm>
#include <utility>

#include <boost/array.hpp>
#include <boost/assert.hpp>

#include "LatticeRealization.hpp"
#include "allowed-momentum.hpp"

template<std::size_t DIM>
class _MomentaSortFunctionObject
{
public:
    bool operator() (const std::pair<lw_vector<int, MAX_DIMENSION>, real_position_t> &v1,
                     const std::pair<lw_vector<int, MAX_DIMENSION>, real_position_t> &v2) const
        {
            BOOST_ASSERT(v1.first.size() == v2.first.size());
            return (v1.second < v2.second);
        }
};

template<std::size_t DIM>
real_position_t _euclidean_norm (const lw_vector<boost::rational<int>, MAX_DIMENSION> coords, const boost::array<boost::array<real_position_t, DIM>, DIM> &primitive_vectors)
{
    // returns distance squared, just like std::norm() does
    BOOST_ASSERT(coords.size() == DIM);
    real_position_t rv = 0;
    // loop over each dimension
    for (unsigned int i = 0; i < DIM; ++i) {
        const real_position_t coords_i_real = boost::rational_cast<real_position_t>(coords[i]);
        real_position_t x = 0;
        // loop over each primitive vector
        for (unsigned int j = 0; j < DIM; ++j) {
            x += coords_i_real * primitive_vectors[j][i];
        }
        rv += x * x;
    }
    return rv;
}

template<std::size_t DIM>
std::vector<lw_vector<int, MAX_DIMENSION> > lowest_momenta (const LatticeRealization<DIM> &lattice, const BoundaryConditions &bcs, unsigned int count)
{
    std::vector<std::pair<lw_vector<int, MAX_DIMENSION>, real_position_t> > pairs;
    //pairs.reserve(number of primitive cells) // (fixme)

    // loop over all momentum sites (yes, it really is this complicated)
    lw_vector<int, MAX_DIMENSION> momentum_site(DIM, 0);
    int not_done = 0;
    boost::array<int*, DIM> overflow_target;
    for (unsigned int i = 0; i < DIM - 1; ++i)
        overflow_target[i] = &momentum_site[i + 1];
    overflow_target[DIM - 1] = &not_done;
    while (not_done == 0) {
#if defined(DEBUG_VMC_LOWEST_MOMENTA) || defined(DEBUG_VMC_ALL)
        for (unsigned int i = 0; i < DIM; ++i)
            std::cerr << momentum_site[i] << ' ';
        std::cerr << std::endl;
#endif

        // actual inner loop (what we're trying to accomplish)
        {
            // determine the norm of the momentum at the site, and save it away (in `pairs`)
            lw_vector<boost::rational<int>, MAX_DIMENSION> momentum(allowed_momentum(momentum_site, lattice, bcs));
            lattice.map_momentum_to_brillouin_zone(momentum);
            const real_position_t euc_norm = _euclidean_norm(momentum, lattice.reciprocal_primitive_vectors);
            pairs.push_back(std::make_pair(momentum_site, euc_norm));
#if defined(DEBUG_VMC_LOWEST_MOMENTA) || defined(DEBUG_VMC_ALL)
            std::cerr << euc_norm << std::endl;
#endif
        }

        // move to next momentum site
        ++momentum_site[0];
        for (unsigned int i = 0; i < DIM; ++i) {
            if (momentum_site[i] == lattice.dimensions[i]) {
                momentum_site[i] = 0;
                ++(*overflow_target[i]);
            } else {
                break;
            }
        }
    }

    // sort things to find the `count` smallest momenta in no particular order
    std::nth_element(pairs.begin(), pairs.begin() + count, pairs.end(), _MomentaSortFunctionObject<DIM>());

    // construct the std::vector for returning
    std::vector<lw_vector<int, MAX_DIMENSION> > rv;
    rv.reserve(count);
    for (unsigned int i = 0; i < count; ++i)
        rv.push_back(pairs[i].first);

#if defined(DEBUG_VMC_LOWEST_MOMENTA) || defined(DEBUG_VMC_ALL)
    std::cerr << "Returning " << rv.size() << " orbitals:" << std::endl;
    for (unsigned int i = 0; i < rv.size(); ++i) {
        for (unsigned int j = 0; j < DIM; ++j)
            std::cerr << rv[i][j] << ' ';
        std::cerr << std::endl;
    }
#endif

    return rv;
}

#endif
