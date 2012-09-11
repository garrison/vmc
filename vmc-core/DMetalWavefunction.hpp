#ifndef _D_METAL_WAVEFUNCTION_HPP
#define _D_METAL_WAVEFUNCTION_HPP

#include <boost/shared_ptr.hpp>

#include "Wavefunction.hpp"
#include "PositionArguments.hpp"

#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * D-wave metal wave function
 *
 * (product of two determinants represents the charge sector with a DBL, then two more
 * determinants for spin up and spin down)
 */
class DMetalWavefunction : public Wavefunction
{
public:
    const boost::shared_ptr<const OrbitalDefinitions> orbital_d1, orbital_d2, orbital_f_up, orbital_f_down;
    const real_t d1_exponent, d2_exponent, f_up_exponent, f_down_exponent;

    DMetalWavefunction (const boost::shared_ptr<const OrbitalDefinitions> &orbital_d1_,
                        const boost::shared_ptr<const OrbitalDefinitions> &orbital_d2_,
                        const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_up_,
                        const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_down_,
                        real_t d1_exponent_,
                        real_t d2_exponent_,
                        real_t f_up_exponent_,
                        real_t f_down_exponent_)
        : Wavefunction(orbital_d1_->get_lattice_ptr()),
          orbital_d1(orbital_d1_),
          orbital_d2(orbital_d2_),
          orbital_f_up(orbital_f_up_),
          orbital_f_down(orbital_f_down_),
          d1_exponent(d1_exponent_),
          d2_exponent(d2_exponent_),
          f_up_exponent(f_up_exponent_),
          f_down_exponent(f_down_exponent_)
        {
        }

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        CeperleyMatrix<amplitude_t> m_cmat_d1, m_cmat_d2, m_cmat_f_up, m_cmat_f_down;
        int m_partial_update_step;

        // the following variables only need be set when a move is in progress
        unsigned int m_up_particles_in_progress, m_down_particles_in_progress;
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const DMetalWavefunction> &wf, const PositionArguments &r_);

    private:
        void perform_move_ (const Move &move);

        template <bool first_pass>
        void do_perform_move (const Move &move);

        amplitude_t psi_ (void) const;

        void finish_move_ (void);

        void cancel_move_ (void);

        void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species);

        void reset_ (const PositionArguments &r_);

        boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const;

        void reset_with_random_configuration (RandomNumberGenerator &rng);

        void reinitialize (void);
    };
};

#endif
