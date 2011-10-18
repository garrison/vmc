#ifndef _SUBSYSTEM_HPP
#define _SUBSYSTEM_HPP

#include <algorithm>

template<class T>
class Subsystem
// this is an abstract base class
{
public:
    virtual bool particle_is_within (const typename T::position_t &position) const = 0;

    virtual ~Subsystem (void)
	{
	}
};

template <class T>
// operates in O(n) time;
// returns a pair, with number of particles in partition being the first thing, parity the second
std::pair<int, int> _fermion_partition (typename T::Arguments &fermions, const Subsystem<T> &subsystem)
{
    // sorts such that everything in the subsystem is put first; everything
    // else after it
    unsigned int i = 0, j = fermions.size() - 1;
    int parity = 1;
    for (;;) {
	if (i >= j)
	    goto done;

	// find a pair of particles that are out of place
	while (subsystem.particle_is_within(fermions[i])) {
	    ++i;
	    if (i == j)
		goto done;
	}
	while (!subsystem.particle_is_within(fermions[j])) {
	    --j;
	    if (i == j)
		// we know that index (i - 1) is the final particle in the subsystem
		return std::pair<int, int>(i, parity);
	}

	// perform swap
	fermions.swap_positions(i, j);
	parity = -parity;

	++i;
	--j;
    }
done:
    // either fermions[i-1] or fermions[i] might be the final particle in the subsystem
    if (i < fermions.size() && subsystem.particle_is_within(fermions[i]))
	++i;
    return std::pair<int, int>(i, parity);
}

template <class T>
std::pair<int, int> fermion_partition (typename T::Arguments &fermions, const Subsystem<T> &subsystem)
{
    std::pair<int, int> rv = _fermion_partition(fermions, subsystem);
#ifdef DEBUG
    unsigned int i = 0;
    for (; i < fermions.size() && subsystem.particle_is_within(fermions[i]); ++i);
    BOOST_ASSERT((int)i == rv.first);
#endif
    return rv;
}


#endif
