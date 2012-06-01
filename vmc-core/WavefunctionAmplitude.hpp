#ifndef _WAVEFUNCTION_AMPLITUDE_HPP
#define _WAVEFUNCTION_AMPLITUDE_HPP

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "PositionArguments.hpp"
#include "Lattice.hpp"
#include "vmc-typedefs.hpp"

class RandomNumberGenerator;

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
     * moves are not allowed until finish_move() has been called.
     */
    void perform_move (Particle particle, unsigned int new_site_index)
        {
            BOOST_ASSERT(!move_in_progress);
            BOOST_ASSERT(r.particle_is_valid(particle));
            BOOST_ASSERT(!r.is_occupied(new_site_index, particle.species) || r[particle] == new_site_index);

            current_move_particle = particle;
            current_move_old_site_index = r[particle];

#if defined(DEBUG_VMC_WAVEFUNCTION_AMPLITUDE) || defined(DEBUG_VMC_ALL)
            if (r[particle] == new_site_index)
                std::cerr << "performing a no-op particle move" << std::endl;
#endif

            perform_move_(particle, new_site_index);

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
    void finish_move (void)
        {
            BOOST_ASSERT(move_in_progress);
            finish_move_();
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = false;
#endif
        }

    /**
     * Cancels the current move, such that new moves are allowed
     */
    void cancel_move (void)
        {
            BOOST_ASSERT(move_in_progress);
            cancel_move_(current_move_particle, current_move_old_site_index);
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = false;
#endif
        }

    void swap_particles (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
        {
            BOOST_ASSERT(!move_in_progress);
            BOOST_ASSERT(species < r.get_N_species());
            BOOST_ASSERT(particle1_index < r.get_N_filled(species));
            BOOST_ASSERT(particle2_index < r.get_N_filled(species));
            BOOST_ASSERT(particle1_index != particle2_index);
            r.swap_particles(particle1_index, particle2_index, species);
            swap_particles_(particle1_index, particle2_index, species);
        }

    /**
     * Resets the wavefunction amplitude with a new set of positions
     */
    void reset (const PositionArguments &r_)
        {
            BOOST_ASSERT(!move_in_progress);
            reset_(r_);
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

    /**
     * Reset the wavefunction's positions randomly.
     *
     * Subclasses may wish to override this method if they have special
     * projection properties, e.g., no spin-up and spin-down particle allowed
     * on the same site.
     */
    virtual void reset_with_random_positions (RandomNumberGenerator &rng);

private:
    /**
     * This method is responsible for updating the PositionArguments "r" as
     * well as updating the state of the object such that psi_() returns the
     * new amplitude.
     */
    virtual void perform_move_ (Particle particle, unsigned int new_site_index) = 0;

    virtual amplitude_t psi_ (void) const = 0;

    virtual void finish_move_ (void) = 0;

    /**
     * This method is responsible for returning the particle (in the
     * PositionArguments "r") to its original position, as well as restoring
     * the state of everything to the way it was before.
     */
    virtual void cancel_move_ (Particle particle, unsigned int old_site_index) = 0;

    // this gets called *after* the particles have been updated in this->r
    virtual void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species) = 0;

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

    // we remember these things for when we want to cancel a move
    Particle current_move_particle;
    unsigned int current_move_old_site_index;

    const boost::shared_ptr<const Lattice> lattice;

private:
#ifndef BOOST_DISABLE_ASSERTS
    bool move_in_progress;
#endif
};

#endif
