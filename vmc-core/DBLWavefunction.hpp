#ifndef _D_B_L_WAVEFUNCTION_HPP
#define _D_B_L_WAVEFUNCTION_HPP

#include <boost/shared_ptr.hpp>

#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * D-wave Bose Liquid wave function
 *
 * (product of two determinants)
 */
class DBLWavefunction : public Wavefunction
{
public:
    const boost::shared_ptr<const OrbitalDefinitions> orbital_def1, orbital_def2;
    const real_t d1_exponent, d2_exponent;

    DBLWavefunction (const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_2, real_t d1_exponent_, real_t d2_exponent_)
        : Wavefunction(orbital_def_1->get_lattice_ptr()),
          orbital_def1(orbital_def_1),
          orbital_def2(orbital_def_2),
          d1_exponent(d1_exponent_),
          d2_exponent(d2_exponent_)
        {
        }

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        CeperleyMatrix<amplitude_t> cmat1, cmat2;

        int m_partial_update_step;

        // the following variable need be set only when a move is in progress
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<DBLWavefunction> &wf_, const PositionArguments &r_);

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

        void reinitialize (void);
    };
};

#endif
