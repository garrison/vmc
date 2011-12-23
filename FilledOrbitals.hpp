#ifndef _FILLED_ORBITALS_HPP
#define _FILLED_ORBITALS_HPP

#include <cstddef>
#include <cmath>
#include <vector>
#include <set>

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/rational.hpp>
#include <boost/math/constants/constants.hpp>

#include "NDLattice.hpp"
#include "OrbitalDefinitions.hpp"
#include "allowed-momentum.hpp"

/**
 * OrbitalDefinitions based on filled momenta and some given boundary conditions
 */
template<std::size_t DIM>
class FilledOrbitals : public OrbitalDefinitions
{
public:
    /**
     * Constructor
     *
     * @param momentum_sites represents which orbitals are filled
     *
     * @param lattice_ the lattice
     *
     * @param bcs the boundary conditions
     */
    FilledOrbitals (const std::vector<boost::array<int, DIM> > &momentum_sites, const boost::shared_ptr<const NDLattice<DIM> > &lattice_, const typename NDLattice<DIM>::BoundaryConditions &bcs)
        : OrbitalDefinitions(momentum_sites.size(), lattice_),
          boundary_conditions(bcs)
        {
            BOOST_ASSERT(momentum_sites.size() <= lattice->total_sites());

#ifndef BOOST_DISABLE_ASSERTS
            // use a set to check that we don't specify a single orbital twice.
            std::set<boost::array<int, DIM> > existing_momentum_sites;
#endif

            for (unsigned int i = 0; i < momentum_sites.size(); ++i) {
                BOOST_ASSERT(existing_momentum_sites.insert(momentum_sites[i]).second); // true if an element was indeed inserted

                boost::array<boost::rational<int>, DIM> momentum = allowed_momentum(momentum_sites[i], *lattice_, bcs);

                // fixme: there ought to be a more efficient way of iterating
                // over lattice sites than using site_from_index() repeatedly
                const complex_t two_pi_i = 2 * boost::math::constants::pi<real_t>() * complex_t(0, 1);
                for (unsigned int j = 0; j < lattice->total_sites(); ++j) {
                    const typename NDLattice<DIM>::Site site = lattice_->site_from_index(j);
                    real_t dot_product = 0;
                    for (unsigned int k = 0; k < DIM; ++k)
                        dot_product += boost::rational_cast<real_t>(momentum[k] * site[k]);
                    orbitals(i, j) = std::exp(two_pi_i * dot_product);
                }
            }
        }

    /**
     * The boundary conditions with which the object was initialized
     */
    const typename NDLattice<DIM>::BoundaryConditions boundary_conditions;
};

#endif
