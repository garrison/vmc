#ifndef _D_METAL_WAVEFUNCTION_AMPLITUDE_HPP
#define _D_METAL_WAVEFUNCTION_AMPLITUDE_HPP

#include <boost/shared_ptr.hpp>

#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"

#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * D-wave metal wave function
 *
 * (product of two determinants represents the charge sector with a DBL, then two more
 * determinants for spin up and spin down)
 */
class DMetalWavefunctionAmplitude : public WavefunctionAmplitude
{
private:
    CeperleyMatrix<amplitude_t> m_cmat_d1, m_cmat_d2, m_cmat_f_up, m_cmat_f_down;
    const boost::shared_ptr<const OrbitalDefinitions> m_orbital_d1, m_orbital_d2, m_orbital_f_up, m_orbital_f_down;
    const real_t m_d1_exponent, m_d2_exponent, m_f_up_exponent, m_f_down_exponent;
    int m_partial_update_step;

    // the following variables only need be set when a move is in progress
    unsigned int m_up_particles_in_progress, m_down_particles_in_progress;
    Move m_current_move;

public:
    DMetalWavefunctionAmplitude (const PositionArguments &r_,
                                 const boost::shared_ptr<const OrbitalDefinitions> &orbital_d1,
                                 const boost::shared_ptr<const OrbitalDefinitions> &orbital_d2,
                                 const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_up,
                                 const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_down,
                                 real_t d1_exponent,
                                 real_t d2_exponent,
                                 real_t f_up_exponent,
                                 real_t f_down_exponent);

private:
    void perform_move_ (const Move &move);

    template <bool first_pass>
    void do_perform_move (const Move &move);

    amplitude_t psi_ (void) const;

    void finish_move_ (void);

    void cancel_move_ (void);

    void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species);

    void reset_ (const PositionArguments &r_);

    boost::shared_ptr<WavefunctionAmplitude> clone_ (void) const;

    void reset_with_random_configuration (RandomNumberGenerator &rng);

    void reinitialize (void);
};

#endif
