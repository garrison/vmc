#ifndef _FREE_FERMION_WAVEFUNCTION_AMPLITUDE_HPP
#define _FREE_FERMION_WAVEFUNCTION_AMPLITUDE_HPP

#include <boost/shared_ptr.hpp>

#include "WavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * Free fermion wave function
 *
 * (a single determinant with no Jastrow factor)
 */
class FreeFermionWavefunctionAmplitude : public WavefunctionAmplitude
{
private:
    CeperleyMatrix<amplitude_t> cmat;
    const boost::shared_ptr<const OrbitalDefinitions> orbital_def;

public:
    FreeFermionWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_);

private:
    void perform_move_ (const Move &move);

    amplitude_t psi_ (void) const;

    void finish_move_ (void);

    void cancel_move_ (void);

    void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species);

    void reset_ (const PositionArguments &r_);

    boost::shared_ptr<WavefunctionAmplitude> clone_ (void) const;

    void reinitialize (void);
};

#endif
