#ifndef _WAVEFUNCTION_HPP
#define _WAVEFUNCTION_HPP

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "PositionArguments.hpp"
#include "Lattice.hpp"
#include "Move.hpp"
#include "vmc-typedefs.hpp"

class RandomNumberGenerator;

/**
 * Abstract base class representing a wavefunction
 *
 * Each Wavefunction object should be invariant after instantiation.
 */
class Wavefunction
{
public:
    const boost::shared_ptr<const Lattice> lattice;

protected:
    Wavefunction (const boost::shared_ptr<const Lattice> &lattice_)
        : lattice(lattice_)
        {
        }

public:
    virtual ~Wavefunction (void)
        {
        }

    /**
     * Abstract base class representing a wavefunction amplitude, which varies
     * as the particles are moved around.
     */
    class Amplitude
    {
    public:
        virtual ~Amplitude (void)
            {
            }

        /**
         * Moves a particle to a new site
         *
         * After this is called, psi() will return the new amplitude, but
         * further moves are not allowed until finish_move() has been called.
         *
         * A null move should never be sent to this method.
         */
        void perform_move (const Move &move);

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
#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
                move_in_progress = false;
#endif
            }

        /**
         * Cancels the current move, such that new moves are allowed
         */
        void cancel_move (void);

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
        boost::shared_ptr<Wavefunction::Amplitude> clone (void) const
            {
                return clone_();
            }

        const Lattice & get_lattice (void) const
            {
                return *wf->lattice;
            }

        /**
         * Returns the current positions of the particles
         */
        const PositionArguments & get_positions (void) const
            {
                return r;
            }

        /**
         * Propose a move.
         *
         * This does not actually perform the move.  Just proposes a sensible
         * one randomly.
         *
         * Must maintain balance.
         *
         * It is allowed for this to return the null move occasionally (which
         * may be useful to obtain balance).
         */
        virtual Move propose_move (RandomNumberGenerator &rng) const;

    private:
        /**
         * This method is responsible for updating the state of the object such
         * that psi_() returns the new amplitude.  The PositionArguments "r"
         * will be updated before this method is called.  Also, a null move
         * should never be sent to this method.
         */
        virtual void perform_move_ (const Move &move) = 0;

        virtual amplitude_t psi_ (void) const = 0;

        virtual void finish_move_ (void) = 0;

        /**
         * This method is responsible for restoring the state of everything to
         * the way it was before.  After this is called, the PositionArguments
         * r will be reverted to its original state.
         */
        virtual void cancel_move_ (void) = 0;

        // this gets called *after* the particles have been updated in this->r
        virtual void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species) = 0;

        virtual void reset_ (const PositionArguments &r_) = 0;

        virtual boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const = 0;

    protected:
        Amplitude (const boost::shared_ptr<const Wavefunction> &wf_, const PositionArguments &r_)
            : wf(wf_),
              r(r_)
#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
        , move_in_progress(false)
#endif
            {
            }

        const Move & get_reverse_move (void) const
            {
                return reverse_move;
            }

        const boost::shared_ptr<const Wavefunction> wf;
        PositionArguments r;

    private:
        // for when we need to cancel a move
        Move reverse_move;

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
        bool move_in_progress;
#endif
    };

    /**
     * Create a corresponding Wavefunction::Amplitude object with nonzero
     * amplitude
     *
     * Subclasses may wish to override this method if they have special
     * projection properties, e.g., no spin-up and spin-down particle allowed
     * on the same site.
     */
    virtual boost::shared_ptr<Wavefunction::Amplitude> create_nonzero_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts=1000000) const;

    virtual boost::shared_ptr<Wavefunction::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, const PositionArguments &r) const = 0;

    virtual unsigned int get_N_species (void) const = 0;

    virtual unsigned int get_N_filled (unsigned int species) const = 0;
};

#endif
