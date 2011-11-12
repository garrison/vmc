#ifndef _SUBSYSTEM_HPP
#define _SUBSYSTEM_HPP

class Lattice;

class Subsystem
// this is an abstract base class
{
public:
    virtual bool particle_is_within (unsigned int site_index, const Lattice &lattice) const = 0;

    virtual bool lattice_makes_sense (const Lattice &lattice) const = 0;

    virtual ~Subsystem (void)
	{
	}
};

static inline int calculate_subsystem_particle_change (const Subsystem &subsystem, unsigned int current_position, unsigned int proposed_position, const Lattice &lattice)
{
    return ((subsystem.particle_is_within(current_position, lattice) ? -1 : 0)
	    + (subsystem.particle_is_within(proposed_position, lattice) ? 1 : 0));
}

#endif
