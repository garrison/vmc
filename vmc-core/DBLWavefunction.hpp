#ifndef _VMC_D_B_L_WAVEFUNCTION_HPP
#define _VMC_D_B_L_WAVEFUNCTION_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"

/**
 * D-wave Bose Liquid wave function
 *
 * (product of two determinants)
 */
template <typename AmplitudeType>
class DBLWavefunction : public Wavefunction<AmplitudeType>
{
public:
    const boost::shared_ptr<const OrbitalDefinitions> orbital_def1, orbital_def2;
    const real_t d1_exponent, d2_exponent;

    DBLWavefunction (const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_2, real_t d1_exponent_, real_t d2_exponent_)
        : Wavefunction<AmplitudeType>(orbital_def_1->get_lattice_ptr()),
          orbital_def1(orbital_def_1),
          orbital_def2(orbital_def_2),
          d1_exponent(d1_exponent_),
          d2_exponent(d2_exponent_)
        {
        }

    class Amplitude : public Wavefunction<AmplitudeType>::Amplitude
    {
    private:
        CeperleyMatrix<AmplitudeType> cmat1, cmat2;

        int m_partial_update_step;

        // the following variable need be set only when a move is in progress
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const DBLWavefunction> &wf_, const PositionArguments &r_);

    private:
        virtual void perform_move_ (const Move &move) override;

        template <bool first_pass>
        void do_perform_move (const Move &move);

        virtual Big<AmplitudeType> psi_ (void) const override;

        virtual void finish_move_ (void) override;

        virtual void cancel_move_ (void) override;

        virtual void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species) override;

        virtual void reset_ (const PositionArguments &r_) override;

        virtual boost::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> clone_ (void) const override;

        void reinitialize (void);

        virtual void check_for_numerical_error (void) const override;
    };

    virtual boost::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction<AmplitudeType> > &this_ptr, const PositionArguments &r) const override
        {
            BOOST_ASSERT(this == this_ptr.get());
            return boost::make_shared<Amplitude>(boost::dynamic_pointer_cast<const DBLWavefunction>(this_ptr), r);
        }

    virtual unsigned int get_N_species (void) const override
        {
            return 1;
        }

    virtual unsigned int get_N_filled (unsigned int species) const override
        {
            BOOST_ASSERT(species == 0);
            return orbital_def1->get_N_filled();
        }
};

#define VMC_SUPPORTED_TYPE(type) extern template class DBLWavefunction<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_TYPE

#endif
