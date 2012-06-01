#ifndef _FILLED_ORBITALS_HPP
#define _FILLED_ORBITALS_HPP

#include <cstddef>
#include <cmath>
#include <vector>
#include <set>

#include <boost/assert.hpp>
#include <boost/rational.hpp>
#include <boost/math/constants/constants.hpp>

#include "Lattice.hpp"
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
    FilledOrbitals (const std::vector<lw_vector<int, MAX_DIMENSION> > &momentum_sites, const boost::shared_ptr<const Lattice> &lattice_, const BoundaryConditions &bcs)
        : OrbitalDefinitions(perform_filling(momentum_sites, *lattice_, bcs), lattice_),
          boundary_conditions(bcs)
        {
        }

    static Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> perform_filling (const std::vector<lw_vector<int, MAX_DIMENSION> > &momentum_sites, const Lattice &lattice, const BoundaryConditions &bcs)
        {
            BOOST_ASSERT(lattice.n_dimensions() == bcs.size());

            const unsigned int N_filled = momentum_sites.size();
            const unsigned int N_sites = lattice.total_sites();

            BOOST_ASSERT(N_filled <= N_sites);

            Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> orbitals(N_filled, N_sites);

            const real_t normalization = std::pow(real_t(N_filled), -.4);

#ifndef BOOST_DISABLE_ASSERTS
            // use a set to check that we don't specify a single orbital twice.
            std::set<lw_vector<int, MAX_DIMENSION> > existing_momentum_sites;
#endif

            for (unsigned int i = 0; i < N_filled; ++i) {
                BOOST_ASSERT(momentum_sites[i].size() == lattice.n_dimensions());
                BOOST_ASSERT(existing_momentum_sites.insert(momentum_sites[i]).second); // true if an element was indeed inserted

                lw_vector<boost::rational<int>, MAX_DIMENSION> momentum = allowed_momentum<DIM>(momentum_sites[i], lattice, bcs);

                // fixme: there ought to be a more efficient way of iterating
                // over lattice sites than using site_from_index() repeatedly
                const complex_t two_pi_i = 2 * boost::math::constants::pi<real_t>() * complex_t(0, 1);
                for (unsigned int j = 0; j < N_sites; ++j) {
                    const LatticeSite site(lattice.site_from_index(j));
                    real_t dot_product = 0;
                    for (unsigned int k = 0; k < lattice.n_dimensions(); ++k)
                        dot_product += boost::rational_cast<real_t>(momentum[k] * site[k]);
                    orbitals(i, j) = std::exp(two_pi_i * dot_product) * normalization;
                }
            }

            return orbitals;
        }

    /**
     * The boundary conditions with which the object was initialized
     */
    const BoundaryConditions boundary_conditions;
};

#endif
