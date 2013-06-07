#ifndef _VMC_LATTICE_HPP
#define _VMC_LATTICE_HPP

#include <cstddef>
#include <utility>
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
    int basis_index = 0;

public:
    // CYTHON-LIMITATION: must have a default constructor, hence default
    // n_dimensions=0 and set_n_dimensions() below
    explicit LatticeSite (unsigned int n_dimensions=0)
        : bs(n_dimensions)
        {
        }

    /**
     * returns the site of the underlying Bravais lattice
     */
    const BravaisSite & bravais_site (void) const
        {
            return bs;
        }

    void set_n_dimensions (unsigned int n_dimensions)
        {
            bs.resize(n_dimensions);
        }

    unsigned int n_dimensions (void) const
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

    // CYTHON-LIMITATION:
    // http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#c-left-values
    void set_bs_coordinate (std::size_t index, int v)
        {
            bs[index] = v;
        }

    bool operator== (const LatticeSite &other) const
        {
            return (bs == other.bs) && (basis_index == other.basis_index);
        }

    bool operator!= (const LatticeSite &other) const
        {
            return (bs != other.bs) || (basis_index != other.basis_index);
        }

    bool operator< (const LatticeSite &other) const;
};

/**
 * N-dimensional lattice
 *
 * This class represents an N-dimensional lattice, where N is specified.  Most
 * of the operations we do don't depend on the particular primitive vectors of
 * the Bravais lattice, so Lattice turns out to be immensely useful, even
 * though it contains no knowledge of the primitive vectors.
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
    Lattice (const lw_vector<int, MAX_DIMENSION> &dimensions_, int basis_indices_=1);

    Lattice (const std::initializer_list<int> &dimensions_, int basis_indices_=1);

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
     * Returns the total number of Bravais sites on the lattice
     */
    unsigned int total_bravais_sites (void) const
        {
            return m_total_sites / basis_indices;
        }

    /**
     * Maps a lattice index (0..N-1) to a LatticeSite (i.e. its coordinates)
     *
     * @see index()
     */
    LatticeSite operator[] (unsigned int n) const;

    /**
     * Maps a LatticeSite to its corresponding lattice index (0..N-1)
     *
     * @see operator[]()
     */
    unsigned int index (const LatticeSite &site) const;

    /**
     * Returns true if the given site is on the lattice
     *
     * Before calling this method, one should ensure that the LatticeSite has
     * the correct number of dimensions such that it would even make sense on
     * this lattice.
     */
    bool site_is_valid (const LatticeSite &site) const;

    /**
     * Adds to a LatticeSite the vector corresponding to the given BravaisSite
     *
     * Addition is done in place.  (The prefix "asm" is meant to remind of
     * this, since addition is typically done in-place in assembly language.)
     *
     * This may return a site that is off the lattice.  You may wish to call
     * enforce_boundary() after calling this.
     *
     * @see asm_subtract_site_vector()
     */
    void asm_add_site_vector (LatticeSite &site, const BravaisSite &other) const;

    /**
     * Subtracts from a LatticeSite the vector corresponding to the given BravaisSite
     *
     * Subtraction is done in place.  (The prefix "asm" is meant to remind of
     * this, since addition is typically done in-place in assembly language.)
     *
     * This may return a site that is off the lattice.  You may wish to call
     * enforce_boundary() after calling this.
     *
     * @see asm_add_site_vector()
     */
    void asm_subtract_site_vector (LatticeSite &site, const BravaisSite &other) const;

    /**
     * If the site is outside the lattice, move it to the corresponding site
     * inside the lattice
     *
     * @return the phase change due to any crossings of the boundary.
     *
     * If the original site is off the lattice in a direction with open
     * boundary conditions, the `site` will still be moved to a site on the
     * lattice, but the phase difference returned will be zero.
     *
     * If the BoundaryConditions array has zero elements, then the boundary
     * will be enforced but it is assumed that the phase is irrelevant.  The
     * value 1 will be returned for the phase in this case, i.e. everything
     * will be treated as if it has periodic boundary conditions.
     */
    template <typename PhaseType>
    PhaseType enforce_boundary (LatticeSite &site, const BoundaryConditions<PhaseType> &bcs) const;

    /**
     * If the site is outside the lattice, move it to the corresponding site
     * inside the lattice assuming periodic boundary conditions
     */
    void enforce_periodic_boundary (LatticeSite &site) const;

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
            enforce_periodic_boundary(site);
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

template <typename PhaseType>
PhaseType Lattice::enforce_boundary (LatticeSite &site, const BoundaryConditions<PhaseType> &bcs) const
{
    BOOST_ASSERT(site.n_dimensions() == n_dimensions());
    BOOST_ASSERT(bcs.size() == 0 || bcs.size() == n_dimensions());
    using std::conj;
    PhaseType phase_change = 1;
    for (unsigned int dim = 0; dim < n_dimensions(); ++dim) {
        while (site[dim] >= dimensions[dim]) {
            site[dim] -= dimensions[dim];
            if (bcs.size() != 0)
                phase_change *= bcs[dim].phase();
        }
        while (site[dim] < 0) {
            site[dim] += dimensions[dim];
            if (bcs.size() != 0)
                phase_change *= conj(bcs[dim].phase());
        }
    }

    // this is often unnecessary ... should it be in a separate
    // function to be called before this one when needed?
    do_safe_modulus(site.basis_index, basis_indices);

    BOOST_ASSERT(site_is_valid(site));
    return phase_change;
}

#endif
