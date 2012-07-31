#ifndef _RVB_WAVEFUNCTION_AMPLITUDE_HPP
#define _RVB_WAVEFUNCTION_AMPLITUDE_HPP

#include <vector>

#include "WavefunctionAmplitude.hpp"
#include "CeperleyMatrix.hpp"

/**
 * RVB (projected BCS wave function at half-filling)
 */
class RVBWavefunctionAmplitude : public WavefunctionAmplitude
{
private:
    // Because we don't yet have a way of doing row-column updates in a single
    // step, on each step we copy m_cmat to m_new_cmat, update it twice, and
    // copy it back if we choose to accept the move.
    bool m_update_in_progress;
    CeperleyMatrix<amplitude_t> m_cmat, m_new_cmat;
    const std::vector<complex_t> m_phi;

public:
    RVBWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const Lattice> &lattice_, const std::vector<complex_t> &phi);

    virtual Move propose_move (RandomNumberGenerator &rng) const;

private:
    void perform_move_ (const Move &move);

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
