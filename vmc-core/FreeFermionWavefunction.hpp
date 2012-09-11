#ifndef _FREE_FERMION_WAVEFUNCTION_HPP
#define _FREE_FERMION_WAVEFUNCTION_HPP

#include <boost/shared_ptr.hpp>

#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * Free fermion wave function
 *
 * (a single determinant with no Jastrow factor)
 */
class FreeFermionWavefunction : public Wavefunction
{
public:
    const boost::shared_ptr<const OrbitalDefinitions> orbital_def;

    FreeFermionWavefunction (const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_)
        : Wavefunction(orbital_def_->get_lattice_ptr()),
          orbital_def(orbital_def_)
        {
        }

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        CeperleyMatrix<amplitude_t> cmat;

    public:
        Amplitude (const boost::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_);

    private:
        void perform_move_ (const Move &move);

        amplitude_t psi_ (void) const;

        void finish_move_ (void);

        void cancel_move_ (void);

        void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species);

        void reset_ (const PositionArguments &r_);

        boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const;

        void reinitialize (void);
    };
};

#endif
