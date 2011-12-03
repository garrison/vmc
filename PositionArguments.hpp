#ifndef _POSITION_ARGUMENTS_HPP
#define _POSITION_ARGUMENTS_HPP

#include <vector>
#include <cstddef>

#include <boost/assert.hpp>

/**
 * Stores each particle's position such that we can also quickly query to see
 * whether a given site is vacant or not.
 *
 * Remember: this object must support having multiple particles on the same
 * site, since this occurs e.g. while SwappedSystem is initializing itself.
 */
class PositionArguments
{
private:
    std::vector<unsigned int> r; /**< stores the position of each particle */
    std::vector<int> positions; /**< stores the particle count on each site */
public:
    PositionArguments (const std::vector<unsigned int> &r_, unsigned int N_sites);

    /**
     * returns the site index of the given particle
     */
    unsigned int operator[] (unsigned int particle) const
        {
            BOOST_ASSERT(particle < get_N_filled());
            return r[particle];
        }

    /**
     * returns the number of particles
     */
    std::size_t size (void) const
        {
            return r.size();
        }

    /** moves a particle to a new position
     *
     * @param particle particle index
     * @param position position index of target site
     */
    void update_position (unsigned int particle, unsigned int position)
        {
            BOOST_ASSERT(particle < get_N_filled());
            BOOST_ASSERT(position < get_N_sites());
            --positions[r[particle]];
            ++positions[position];
            r[particle] = position;
        }

    bool is_occupied (unsigned int position) const
        {
            return positions[position] != 0;
        }

    unsigned int get_N_sites (void) const
        {
            return positions.size();
        }

    unsigned int get_N_filled (void) const
        {
            return r.size();
        }
};

#endif
