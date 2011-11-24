#ifndef _WAVEFUNCTION_AMPLITUDE_HPP
#define _WAVEFUNCTION_AMPLITUDE_HPP

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "PositionArguments.hpp"
#include "Lattice.hpp"
#include "vmc-typedefs.hpp"

class WavefunctionAmplitude
{
public:
    virtual ~WavefunctionAmplitude (void)
        {
        }

    void move_particle (unsigned int particle, unsigned int new_site_index)
        {
            BOOST_ASSERT(!move_in_progress);
            BOOST_ASSERT(particle < r.get_N_filled());
            BOOST_ASSERT(!r.is_occupied(new_site_index) || r[particle] == new_site_index);
#ifdef DEBUG
            if (r[particle] == new_site_index)
                std::cerr << "performing a no-op particle move" << std::endl;
#endif
            move_particle_(particle, new_site_index);
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = true;
#endif
        }

    amplitude_t psi (void) const
        {
            return psi_();
        }

    void finish_particle_moved_update (void)
        {
            BOOST_ASSERT(move_in_progress);
            finish_particle_moved_update_();
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = false;
#endif
        }

    void reset (const PositionArguments &r_)
        {
            reset_(r_);
#ifndef BOOST_DISABLE_ASSERTS
            move_in_progress = false;
#endif
        }

    boost::shared_ptr<WavefunctionAmplitude> clone (void) const
        {
            return clone_();
        }

    const Lattice & get_lattice (void) const
        {
            return *lattice;
        }

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
