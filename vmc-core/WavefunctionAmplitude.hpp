#ifndef _WAVEFUNCTION_AMPLITUDE_HPP
#define _WAVEFUNCTION_AMPLITUDE_HPP

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "PositionArguments.hpp"
#include "Lattice.hpp"
#include "vmc-typedefs.hpp"

/**
 * Abstract base class representing a wavefunction amplitude
 *
 * Specifically, a WavefunctionAmplitude keeps track of both the wavefunction
 * itself (which is invariant after instantiation) and its amplitude, which
 * varies as the particles are moved around.
 */
class WavefunctionAmplitude
{
public:
    virtual ~WavefunctionAmplitude (void)
        {
        }

    /**
     * Moves a particle to a new site
     *
     * After this is called, psi() will return the new amplitude, but further
     * moves are not allowed until finish_particle_moved_update() has been
     * called.
     */
    void move_particle (unsigned int particle, unsigned int new_site_index)
        {
            BOOST_ASSERT(!move_in_progress);
            BOOST_ASSERT(particle < r.get_N_filled());
            BOOST_ASSERT(!r.is_occupied(new_site_index) || r[particle] == new_site_index);
#if defined(DEBUG_VMC_WAVEFUNCTION_AMPLITUDE) || defined(DEBUG_VMC_ALL)
            if (r[particle] == new_site_index)
                std::cerr << "performing a no-op particle move" << std::endl;
#endif
            move_particle_(particle, new_site_index);
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = true;
#endif
        }

    /**
     * Returns the current amplitude of the wavefunction
     */
    amplitude_t psi (void) const
        {
            return psi_();
        }

    /**
     * Completes the current move, such that new moves are allowed
     */
    void finish_particle_moved_update (void)
        {
            BOOST_ASSERT(move_in_progress);
            finish_particle_moved_update_();
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = false;
#endif
        }

    /**
     * Resets the wavefunction amplitude with a new set of positions
     */
    void reset (const PositionArguments &r_)
        {
            reset_(r_);
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = false;
#endif
        }

    /**
     * Returns a clone of the current object
     */
    boost::shared_ptr<WavefunctionAmplitude> clone (void) const
        {
            return clone_();
        }

    const Lattice & get_lattice (void) const
        {
            return *lattice;
        }

    /**
     * Returns the current positions of the particles
     */
    const PositionArguments & get_positions (void) const
        {
            return r;
        }

private:
    virtual void move_particle_ (unsigned int particle, unsigned int new_site_index) = 0;

    virtual amplitude_t psi_ (void) const = 0;

    virtual void finish_particle_moved_update_ (void) = 0;

    virtual void reset_ (const PositionArguments &r_) = 0;

    virtual boost::shared_ptr<WavefunctionAmplitude> clone_ (void) const = 0;

protected:
    WavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const Lattice> &lattice_)
    : r(r_),
      lattice(lattice_)
#ifndef BOOST_DISABLE_ASSERTS
    , move_in_progress(false)
#endif
        {
        }

    PositionArguments r;

    const boost::shared_ptr<const Lattice> lattice;

private:
#ifndef BOOST_DISABLE_ASSERTS
    bool move_in_progress;
#endif
};

#endif
