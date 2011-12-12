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
    bool m_particle_moved_is_up;

public:
    DMetalWavefunctionAmplitude (const PositionArguments &r_,
                              const boost::shared_ptr<const OrbitalDefinitions> &orbital_d1,
                              const boost::shared_ptr<const OrbitalDefinitions> &orbital_d2,
                              const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_up,
                              const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_down);

private:
    void move_particle_ (unsigned int particle, unsigned int new_site_index);

    amplitude_t psi_ (void) const;

    void finish_particle_moved_update_ (void);

    void reset_ (const PositionArguments &r_);

    boost::shared_ptr<WavefunctionAmplitude> clone_ (void) const;

    void reinitialize (void);
};

#endif
