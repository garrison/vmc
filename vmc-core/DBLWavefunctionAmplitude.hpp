#ifndef _D_B_L_WAVEFUNCTION_AMPLITUDE_HPP
#define _D_B_L_WAVEFUNCTION_AMPLITUDE_HPP

#include <boost/shared_ptr.hpp>

#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * D-wave Bose Liquid wave function
 *
 * (product of two determinants)
 */
class DBLWavefunctionAmplitude : public WavefunctionAmplitude
{
private:
    CeperleyMatrix<amplitude_t> cmat1, cmat2;
    const boost::shared_ptr<const OrbitalDefinitions> orbital_def1, orbital_def2;
    const real_t d1_exponent, d2_exponent;

public:
    DBLWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_2, real_t d1_exponent_, real_t d2_exponent_);

private:
    void move_particle_ (Particle particle, unsigned int new_site_index);

    amplitude_t psi_ (void) const;

    void finish_particle_moved_update_ (void);

    void reset_ (const PositionArguments &r_);

    boost::shared_ptr<WavefunctionAmplitude> clone_ (void) const;

    void reinitialize (void);
};

#endif
