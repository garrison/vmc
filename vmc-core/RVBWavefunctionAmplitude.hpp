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
    CeperleyMatrix<amplitude_t> m_cmat;
    const std::vector<complex_t> m_phi;

public:
    RVBWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const Lattice> &lattice_, const std::vector<complex_t> &phi);

private:
    void perform_move_ (Particle particle, unsigned int new_site_index);

    amplitude_t psi_ (void) const;

    void finish_move_ (void);

    void reset_ (const PositionArguments &r_);

    boost::shared_ptr<WavefunctionAmplitude> clone_ (void) const;

    void reset_with_filler (const RandomFiller &filler, rng_class &rng);

    void reinitialize (void);
};

#endif
