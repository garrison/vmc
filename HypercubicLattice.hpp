#ifndef _HYPERCUBIC_LATTICE_HPP
#define _HYPERCUBIC_LATTICE_HPP

#include <cstddef>

#include <boost/array.hpp>

#include "Lattice.hpp"
#include "BoundaryCondition.hpp"
#include "random-move.hpp"

template<std::size_t DIM>
class HypercubicLattice : public Lattice
{
public:
    typedef boost::array<int, DIM> Site;
    typedef boost::array<BoundaryCondition, DIM> BoundaryConditions;

    static const unsigned int dimensions = DIM;
    static const unsigned int move_axes = DIM;
    static const unsigned int basis_indices = 1; // Bravais lattice

    HypercubicLattice (const boost::array<int, DIM> &length_)
        : Lattice(count_total_sites(length_)),
          length(length_)
        {
        }

    Site site_from_index (unsigned int n) const
        {
            BOOST_ASSERT(n < total_sites());
            Site rv;
            for (unsigned int i = 0; i < DIM; ++i) {
                rv[i] = n % length[i];
                n /= length[i];
            }
            BOOST_ASSERT(site_is_valid(rv));
            return rv;
        }

    unsigned int site_to_index (const Site &site) const
        {
            BOOST_ASSERT(site_is_valid(site));

            unsigned int c = 1, n = 0;
            for (unsigned int i = 0; i < DIM; ++i) {
                n += site[i] * c;
                c *= length[i];
            }
            BOOST_ASSERT(site == site_from_index(n));
            return n;
        }

    bool site_is_valid (const Site &site) const
        {
            for (unsigned int i = 0; i < DIM; ++i) {
                if (site[i] >= length[i])
                    return false;
            }
            return true;
        }

    unsigned int basis_index (const Site &site) const
        {
            // a Bravais lattice always returns zero
            (void) site;
            return 0;
        }

    void move_site (Site &site, unsigned int move_axis, int step_direction) const
        {
            BOOST_ASSERT(move_axis < move_axes);
            BOOST_ASSERT(step_direction == -1 || step_direction == 1);
            site[move_axis] += step_direction;
            // enforce PBC
            if (site[move_axis] >= length[move_axis])
                site[move_axis] -= length[move_axis];
            else if (site[move_axis] < 0)
                site[move_axis] += length[move_axis];
            BOOST_ASSERT(site_is_valid(site));
        }

    phase_t asm_add_site_vector (Site &site, const Site &other, const BoundaryConditions *bcs=0) const
        {
            BOOST_ASSERT(basis_index(other) == 0);
            for (unsigned int i = 0; i < DIM; ++i)
                site[i] += other[i];
            return enforce_boundary(site, bcs);
        }

    phase_t asm_subtract_site_vector (Site &site, const Site &other, const BoundaryConditions *bcs=0) const
        {
            BOOST_ASSERT(basis_index(other) == 0);
            for (unsigned int i = 0; i < DIM; ++i)
                site[i] -= other[i];
            return enforce_boundary(site, bcs);
        }

    phase_t enforce_boundary (Site &site, const BoundaryConditions *bcs=0) const
        {
            phase_t phase_change = 1;
            for (unsigned int dim = 0; dim < DIM; ++dim) {
                while (site[dim] >= length[dim]) {
                    site[dim] -= length[dim];
                    if (bcs)
                        phase_change *= (*bcs)[dim].phase();
                }
                while (site[dim] < 0) {
                    site[dim] += length[dim];
                    if (bcs)
                        phase_change /= (*bcs)[dim].phase();
                }
            }

            BOOST_ASSERT(site_is_valid(site));
            return phase_change;
        }

    Site move_to_basis_index (const Site &site, unsigned int index) const
        {
            // provides access to any other sites within a given unit cell
            BOOST_ASSERT(index < basis_indices);
            // since this is a Bravais lattice, just return the current site
            (void) index;
            return site;
        }

    unsigned int plan_particle_move_to_nearby_empty_site_virtual (unsigned int particle, const PositionArguments &r, rng_class &rng) const
        {
            return plan_particle_move_to_nearby_empty_site_template(particle, r, *this, rng);
        }

private:
    static inline unsigned int count_total_sites (const boost::array<int, DIM> &length)
        {
            unsigned int rv = 1;
            for (unsigned int i = 0; i < DIM; ++i) {
                BOOST_ASSERT(length[i] > 0);
                rv *= length[i];
            }
            return rv;
        }

    const boost::array<int, DIM> length;
};

#endif
