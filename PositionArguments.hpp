#ifndef _POSITION_ARGUMENTS_HPP
#define _POSITION_ARGUMENTS_HPP

#include <vector>
#include <cstddef>

#include <boost/assert.hpp>

// remember: for a brief period of time, this object may contain multiple
// particles on the same site, e.g. while SwappedSystem is initializing itself.

class PositionArguments
// stores each particle's position such that we can also quickly query to see
// whether there is a particle on a given site.
{
private:
    std::vector<unsigned int> r; // stores the position of each particle
    std::vector<int> positions; // stores the particle count on each site
public:
    PositionArguments (const std::vector<unsigned int> &r_, unsigned int N_sites);

    unsigned int operator[] (unsigned int particle) const
	{
	    BOOST_ASSERT(particle < get_N_filled());
	    return r[particle];
	}

    std::size_t size (void) const
	{
	    return r.size();
	}

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
