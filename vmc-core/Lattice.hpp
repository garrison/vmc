#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include "PositionArguments.hpp"
#include "vmc-typedefs.hpp"

/**
 * Abstract base class for a lattice of any dimension.
 *
 * Any code that does not depend on the details of the lattice can use this
 * interface directly.
 */
class Lattice
{
public:
    virtual ~Lattice (void)
        {
        }

    /**
     * Returns the total number of sites on the lattice
     */
    unsigned int total_sites (void) const
        {
            return m_total_sites;
        }

    /**
     * Virtual function, called by plan_particle_move_to_nearby_empty_site()
     *
     * (plan_particle_move_to_nearby_empty_site() in random-move.hpp provides a
     * more natural API, but since the implementation depends on the details of
     * the underlying lattice, our only choice is to implement it as a virtual
     * function, which we do here.)
     *
     * @see plan_particle_move_to_nearby_empty_site()
     */
    virtual unsigned int plan_particle_move_to_nearby_empty_site_virtual (Particle particle, const PositionArguments &r, rng_class &rng) const = 0;

protected:
    Lattice (unsigned int total_sites)
        : m_total_sites(total_sites)
        {
        }

    const unsigned int m_total_sites;
};

#include <cstddef>
#include <vector>

#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random.hpp>

#include "BoundaryCondition.hpp"
#include "vmc-math-utils.hpp"

/**
 * N-dimensional lattice
 *
 * This class represents an N-dimensional lattice, where N is specified.  Most
 * of the operations we do don't depend on the particular primitive vectors of
 * the Bravais lattice, so NDLattice turns out to be immensely useful, even
 * though it contains no knowledge of the primitive vectors.  For physical
 * realizations of lattices in real space, see LatticeRealization and its
 * subclasses, e.g. HypercubicLattice.
 *
 * @see LatticeRealization
 */
template<std::size_t DIM>
class NDLattice : public Lattice
{
public:
    BOOST_STATIC_ASSERT(DIM > 0);

    /** number of dimensions of the lattice */
    static const unsigned int dimensions = DIM;

    /** boundary conditions in each direction */
    typedef boost::array<BoundaryCondition, DIM> BoundaryConditions;

    /**
     * coordinates of a site on the Bravais lattice, represented as an array of
     * zero-based indices
     *
     * The real-space position of a site would be given by the sum of each
     * index multiplied by its respective primitive vector.
     */
    typedef boost::array<int, DIM> BravaisSite;

    /**
     * a site on the lattice, represented by its BravaisSite and basis index
     * (which is always zero if the physical system's sites all live on a
     * Bravais lattice)
     *
     * @see BravaisSite
     */
    struct Site
    {
    private:
        BravaisSite bs;
    public:
        int basis_index;

        /**
         * returns the site of the underlying Bravais lattice
         */
        const BravaisSite & bravais_site (void) const
            {
                return bs;
            }

        int operator[] (std::size_t index) const
            {
                return bs[index];
            }

        int & operator[] (std::size_t index)
            {
                return bs[index];
            }

        bool operator== (const Site &other) const
            {
                return (bs == other.bs) && (basis_index == other.basis_index);
            }

        bool operator!= (const Site &other) const
            {
                return (bs != other.bs) || (basis_index != other.basis_index);
            }
    };

protected:
    /**
     * An axis by which we might want to move in configuration space.
     *
     * Each member here represents some step size.
     */
    struct Move {
        boost::array<int, DIM> bravais_site;
        int basis_index;
    };

public:
    /** NDLattice constructor
     *
     * @param length_ an array representing the number of sites in each dimension
     * @param basis_indices_ the number of sites per unit cell of the Bravais lattice
     */
    NDLattice (const boost::array<int, DIM> &length_, int basis_indices_=1)
        : Lattice(count_total_sites(length_, basis_indices_)),
          length(length_),
          basis_indices(basis_indices_)
        {
            // set up offsets
            unsigned int c = 1;
            for (unsigned int i = 0; i < DIM; ++i) {
                offset[i] = c;
                c *= length[i];
            }
            basis_offset = c;

            // set up default move axes
            for (unsigned int i = 0; i < DIM; ++i) {
                Move m;
                m.bravais_site.assign(0);
                m.bravais_site[i] = 1;
                m.basis_index = 0;
                move_axes.push_back(m);
            }
            if (basis_indices > 1) {
                Move m;
                m.bravais_site.assign(0);
                m.basis_index = 1;
                move_axes.push_back(m);
            }
        }

    /**
     * Maps a lattice index (0..N-1) to a Site (i.e. its coordinates)
     *
     * @see site_to_index()
     */
    Site site_from_index (unsigned int n) const
        {
            BOOST_ASSERT(n < total_sites());
            Site rv;
            for (unsigned int i = 0; i < DIM; ++i) {
                rv[i] = n % length[i];
                n /= length[i];
            }
            rv.basis_index = n;
            BOOST_ASSERT(site_is_valid(rv));
            return rv;
        }

    /**
     * Maps a Site to its corresponding lattice index (0..N-1)
     *
     * @see site_from_index()
     */
    unsigned int site_to_index (const Site &site) const
        {
            BOOST_ASSERT(site_is_valid(site));

            unsigned int n = 0;
            for (unsigned int i = 0; i < DIM; ++i) {
                n += site[i] * offset[i];
            }
            n += site.basis_index * basis_offset;
            BOOST_ASSERT(site == site_from_index(n));
            return n;
        }

    bool site_is_valid (const Site &site) const
        {
            for (unsigned int i = 0; i < DIM; ++i) {
                if (site[i] >= length[i] || site[i] < 0)
                    return false;
            }
            if (site.basis_index >= basis_indices || site.basis_index < 0)
                return false;
            return true;
        }

    /**
     * Adds to a Site the vector corresponding to the given BravaisSite
     *
     * Addition is done in place.  (The prefix "asm" is meant to remind of
     * this, since addition is typically done in-place in assembly language.)
     *
     * @return phase change due to any crossings in the boundary.  If boundary
     * conditions are not given, the phase returned will always be 1.
     *
     * @see asm_subtract_site_vector()
     */
    phase_t asm_add_site_vector (Site &site, const BravaisSite &other, const BoundaryConditions *bcs=0) const
        {
            for (unsigned int i = 0; i < DIM; ++i)
                site[i] += other[i];
            return enforce_boundary(site, bcs);
        }

    /**
     * Subtracts from a Site the vector corresponding to the given BravaisSite
     *
     * Subtraction is done in place.  (The prefix "asm" is meant to remind of
     * this, since addition is typically done in-place in assembly language.)
     *
     * @return phase change due to any crossings in the boundary.  If boundary
     * conditions are not given, the phase returned will always be 1.
     *
     * @see asm_add_site_vector()
     */
    phase_t asm_subtract_site_vector (Site &site, const BravaisSite &other, const BoundaryConditions *bcs=0) const
        {
            for (unsigned int i = 0; i < DIM; ++i)
                site[i] -= other[i];
            return enforce_boundary(site, bcs);
        }

    /**
     * If the site is outside the lattice, move it to the corresponding site
     * inside the lattice
     *
     * @return the phase change due to any crossings of the boundary.  If
     * boundary conditions are not given, the value 1 will be returned.
     */
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

            // this is often unnecessary ... should it be in a separate
            // function to be called before this one when needed?
            do_safe_modulus(site.basis_index, basis_indices);

            BOOST_ASSERT(site_is_valid(site));
            return phase_change;
        }

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
     * Any subclasses of NDLattice (e.g. specific lattice realizations) can add
     * additional move axes during object construction.  It may be possible to
     * get better statistics this way.  (For instance, on a triangular lattice
     * we may wish to attempt moves in three possible directions, even though
     * there are only two primitive vectors.)
     *
     * @see plan_particle_move_to_nearby_empty_site()
     * @see plan_particle_move_to_nearby_empty_site_virtual()
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
    void move_site (typename NDLattice<DIM>::Site &site, unsigned int move_axis, int step_direction) const
        {
            BOOST_ASSERT(move_axis < move_axes.size());
            BOOST_ASSERT(step_direction == -1 || step_direction == 1);
            const Move &m = move_axes[move_axis];
            for (unsigned int i = 0; i < DIM; ++i)
                site[i] += step_direction * m.bravais_site[i];
            site.basis_index += step_direction * m.basis_index;
            enforce_boundary(site);
        }

    /**
     * Returns an empty site index that is "near" a given particle.
     *
     * Typically you should call the corresponding function in random-move.hpp
     * instead, which is a wrapper for this.
     *
     * @see plan_particle_move_to_nearby_empty_site()
     */
    unsigned int plan_particle_move_to_nearby_empty_site_virtual (Particle particle, const PositionArguments &r, rng_class &rng) const
        {
            BOOST_ASSERT(r.particle_is_valid(particle));

            unsigned int move_axis;
            if (this->move_axes_count() == 1) {
                move_axis = 0;
            } else {
                boost::uniform_smallint<> axis_distribution(0, this->move_axes_count() - 1);
                boost::variate_generator<rng_class&, boost::uniform_smallint<> > axis_gen(rng, axis_distribution);
                move_axis = axis_gen();
            }

            boost::uniform_smallint<> direction_distribution(0, 1);
            boost::variate_generator<rng_class&, boost::uniform_smallint<> > direction_gen(rng, direction_distribution);
            int step_direction = direction_gen() * 2 - 1;

            const unsigned int original_site_index = r[particle];
            Site site = this->site_from_index(original_site_index);
            unsigned int site_index;
            do {
                this->move_site(site, move_axis, step_direction);
                site_index = this->site_to_index(site);
            } while (r.is_occupied(site_index, particle.species) && site_index != original_site_index);

            return site_index;
        }

private:
    static inline unsigned int count_total_sites (const boost::array<int, DIM> &length, int basis_indices)
        {
            unsigned int rv = 1;
            for (unsigned int i = 0; i < DIM; ++i) {
                BOOST_ASSERT(length[i] > 0);
                rv *= length[i];
            }
            BOOST_ASSERT(basis_indices > 0);
            rv *= basis_indices;
            return rv;
        }

public:
    const boost::array<int, DIM> length;
    const int basis_indices; /**< number of sites per Bravais unit cell */

private:
    // these both remain constant after initialization as well
    boost::array<int, DIM> offset;
    int basis_offset;

protected:
    // this can be modified at will until the object is fully instantiated, but
    // after that it should not be changed
    std::vector<struct Move> move_axes;
};

#endif
