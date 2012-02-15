#ifndef _SUBSYSTEM_HPP
#define _SUBSYSTEM_HPP

class Lattice;

/**
 * Abstract base class for representing some subset of a lattice's sites
 */
class Subsystem
{
public:
    /**
     * Returns true if the given site index is within the subsystem
     */
    virtual bool position_is_within (unsigned int site_index, const Lattice &lattice) const = 0;

    /**
     * Literally, returns true if the lattice "makes sense" for the current subsystem.
     *
     * The number of dimensions of the lattice needs to be the same as for the
     * defined subsystem, and the subsystem has to "fit" on the given lattice
     * for this to return true.
     */
    virtual bool lattice_makes_sense (const Lattice &lattice) const = 0;

    virtual ~Subsystem (void)
        {
        }
};

class WavefunctionAmplitude;

/**
 * Returns the number of particles in the given subsystem
 */
extern unsigned int count_N_subsystem (const WavefunctionAmplitude &wf, const Subsystem &subsystem);

/**
 * Returns the change in the number of particles in the subsystem if we were to
 * move a particle from current_position to proposed_position.
 *
 * @return 1, 0, or -1
 */
static inline int calculate_subsystem_particle_change (const Subsystem &subsystem, unsigned int current_position, unsigned int proposed_position, const Lattice &lattice)
{
    return ((subsystem.position_is_within(current_position, lattice) ? -1 : 0)
            + (subsystem.position_is_within(proposed_position, lattice) ? 1 : 0));
}

#endif
