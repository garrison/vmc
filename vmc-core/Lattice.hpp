#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include <cstddef>
#include <vector>

#include <boost/assert.hpp>

#include "PositionArguments.hpp"
#include "vmc-typedefs.hpp"
#include "lw_vector.hpp"
#include "BoundaryCondition.hpp"
#include "vmc-math-utils.hpp"

/**
 * coordinates of a site on the Bravais lattice, represented as an array of
 * zero-based indices
 *
 * The real-space position of a site would be given by the sum of each
 * index multiplied by its respective primitive vector.
 */
typedef lw_vector<int, MAX_DIMENSION> BravaisSite;

class Lattice;

/**
 * a site on the lattice, represented by its BravaisSite and basis index
 * (which is always zero if the physical system's sites all live on a
 * Bravais lattice)
 *
 * @see BravaisSite
 */
struct LatticeSite
{
private:
    BravaisSite bs;
public:
    explicit LatticeSite (const Lattice &lattice);

    int basis_index;

    /**
     * returns the site of the underlying Bravais lattice
     */
    const BravaisSite & bravais_site (void) const
        {
            return bs;
        }

    unsigned int size (void) const
        {
            return bs.size();
        }

    int operator[] (std::size_t index) const
        {
            return bs[index];
        }

    int & operator[] (std::size_t index)
        {
            return bs[index];
        }

    bool operator== (const LatticeSite &other) const
        {
            return (bs == other.bs) && (basis_index == other.basis_index);
        }

    bool operator!= (const LatticeSite &other) const
        {
            return (bs != other.bs) || (basis_index != other.basis_index);
        }
};

/**
 * N-dimensional lattice
 *
 * This class represents an N-dimensional lattice, where N is specified.  Most
 * of the operations we do don't depend on the particular primitive vectors of
 * the Bravais lattice, so Lattice turns out to be immensely useful, even
 * though it contains no knowledge of the primitive vectors.  For physical
 * realizations of lattices in real space, see LatticeRealization and its
 * subclasses, e.g. HypercubicLattice.
 *
 * @see LatticeRealization
 */
class Lattice
{
protected:
    /**
     * An axis by which we might want to move in configuration space.
     *
     * Each member here represents some step size.
     */
    struct MoveAxis {
        BravaisSite bravais_site;
        int basis_index;
    };

public:
    /** Lattice constructor
     *
     * @param dimensions_ an array representing the number of sites in each dimension
     * @param basis_indices_ the number of sites per unit cell of the Bravais lattice
     */
    Lattice (const lw_vector<int, MAX_DIMENSION> &dimensions_, int basis_indices_=1)
        : dimensions(dimensions_),
          basis_indices(basis_indices_),
          m_total_sites(count_total_sites(dimensions_, basis_indices_)),
          offset(dimensions_.size())
        {
            BOOST_ASSERT(dimensions.size() > 0);
            BOOST_ASSERT(basis_indices > 0);

            // set up offsets
            unsigned int c = 1;
            for (unsigned int i = 0; i < dimensions.size(); ++i) {
                offset[i] = c;
                c *= dimensions[i];
            }
            basis_offset = c;

            // set up default move axes
            for (unsigned int i = 0; i < dimensions.size(); ++i) {
                MoveAxis m;
                m.bravais_site.resize(dimensions.size(), 0);
                m.bravais_site[i] = 1;
                m.basis_index = 0;
                move_axes.push_back(m);
            }
            if (basis_indices > 1) {
                MoveAxis m;
                m.bravais_site.resize(dimensions.size(), 0);
                m.basis_index = 1;
                move_axes.push_back(m);
            }
        }

    /**
     * Returns the number of dimensions of the lattice
     */
    unsigned int n_dimensions (void) const
        {
            return dimensions.size();
        }

    /**
     * Returns the total number of sites on the lattice
     */
    unsigned int total_sites (void) const
        {
            return m_total_sites;
        }

    /**
     * Maps a lattice index (0..N-1) to a LatticeSite (i.e. its coordinates)
     *
     * @see site_to_index()
     */
    LatticeSite site_from_index (unsigned int n) const;

    /**
     * Maps a LatticeSite to its corresponding lattice index (0..N-1)
     *
     * @see site_from_index()
     */
    unsigned int site_to_index (const LatticeSite &site) const;

    /**
     * Returns true if the given site is on the lattice
     */
    bool site_is_valid (const LatticeSite &site) const;

    /**
     * Adds to a LatticeSite the vector corresponding to the given BravaisSite
     *
     * Addition is done in place.  (The prefix "asm" is meant to remind of
     * this, since addition is typically done in-place in assembly language.)
     *
     * @return phase change due to any crossings in the boundary.  If boundary
     * conditions are not given, the phase returned will always be 1.
     *
     * @see asm_subtract_site_vector()
     */
    phase_t asm_add_site_vector (LatticeSite &site, const BravaisSite &other, const BoundaryConditions *bcs=0) const;

    /**
     * Subtracts from a LatticeSite the vector corresponding to the given BravaisSite
     *
     * Subtraction is done in place.  (The prefix "asm" is meant to remind of
     * this, since addition is typically done in-place in assembly language.)
     *
     * @return phase change due to any crossings in the boundary.  If boundary
     * conditions are not given, the phase returned will always be 1.
     *
     * @see asm_add_site_vector()
     */
    phase_t asm_subtract_site_vector (LatticeSite &site, const BravaisSite &other, const BoundaryConditions *bcs=0) const;

    /**
     * If the site is outside the lattice, move it to the corresponding site
     * inside the lattice
     *
     * @return the phase change due to any crossings of the boundary.  If
     * boundary conditions are not given, the value 1 will be returned.
     */
    phase_t enforce_boundary (LatticeSite &site, const BoundaryConditions *bcs=0) const;

    /**
     * Returns the number of move axes
     *
     * "Move axis" here means a direction in which one can take one or more
     * steps in configuration space.  Before calling
     * plan_particle_move_to_nearby_empty_site(), one first chooses a move
     * axis, typically at random.
     *
     * By default, there is one move axis per dimension in the lattice, plus
     * (if the lattice is not a Bravais lattice) a move axis that only changes
     * the basis_index (and thus moves particles around within a unit cell of
     * the underlying Bravais lattice).
     *
     * Any subclasses of Lattice (e.g. specific lattice realizations) can add
     * additional move axes during object construction.  It may be possible to
     * get better statistics this way.  (For instance, on a triangular lattice
     * we may wish to attempt moves in three possible directions, even though
     * there are only two primitive vectors.)
     *
     * @see plan_particle_move_to_nearby_empty_site()
     */
    unsigned int move_axes_count (void) const
        {
            return move_axes.size();
        }

    /**
     * Moves site a single step in a given direction
     *
     * @param site site to move
     * @param move_axis index representing the axis to move along
     * @param step_direction +1 or -1 depending on which direction to move
     */
    void move_site (LatticeSite &site, unsigned int move_axis, int step_direction) const
        {
            BOOST_ASSERT(site_is_valid(site));
            BOOST_ASSERT(move_axis < move_axes.size());
            BOOST_ASSERT(step_direction == -1 || step_direction == 1);
            const MoveAxis &m = move_axes[move_axis];
            BOOST_ASSERT(m.bravais_site.size() == n_dimensions());
            for (unsigned int i = 0; i < n_dimensions(); ++i)
                site[i] += step_direction * m.bravais_site[i];
            site.basis_index += step_direction * m.basis_index;
            enforce_boundary(site);
        }

private:
    static inline unsigned int count_total_sites (const lw_vector<int, MAX_DIMENSION> &dimensions, int basis_indices)
        {
            BOOST_ASSERT(dimensions.size() > 0);
            unsigned int rv = 1;
            for (unsigned int i = 0; i < dimensions.size(); ++i) {
                BOOST_ASSERT(dimensions[i] > 0);
                rv *= dimensions[i];
            }
            BOOST_ASSERT(basis_indices > 0);
            rv *= basis_indices;
            return rv;
        }

public:
    const lw_vector<int, MAX_DIMENSION> dimensions;
    const int basis_indices; /**< number of sites per Bravais unit cell */

private:
    // these all remain constant after initialization as well
    const unsigned int m_total_sites;
    lw_vector<int, MAX_DIMENSION> offset;
    int basis_offset;

protected:
    // this can be modified at will until the object is fully instantiated, but
    // after that it should not be changed
    std::vector<struct MoveAxis> move_axes;
};

inline LatticeSite::LatticeSite (const Lattice &lattice)
  : bs(lattice.n_dimensions())
{
}

#endif
