#ifndef _VMC_D_METAL_WAVEFUNCTION_HPP
#define _VMC_D_METAL_WAVEFUNCTION_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

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
template <typename AmplitudeType>
class DMetalWavefunction : public Wavefunction<AmplitudeType>
{
public:
    const boost::shared_ptr<const OrbitalDefinitions> orbital_d1, orbital_d2, orbital_f_up, orbital_f_dn;
    const real_t d1_exponent, d2_exponent, f_up_exponent, f_dn_exponent;

    DMetalWavefunction (const boost::shared_ptr<const OrbitalDefinitions> &orbital_d1_,
                        const boost::shared_ptr<const OrbitalDefinitions> &orbital_d2_,
                        const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_up_,
                        const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_dn_,
                        real_t d1_exponent_,
                        real_t d2_exponent_,
                        real_t f_up_exponent_,
                        real_t f_dn_exponent_)
        : Wavefunction<AmplitudeType>(orbital_d1_->get_lattice_ptr()),
          orbital_d1(orbital_d1_),
          orbital_d2(orbital_d2_),
          orbital_f_up(orbital_f_up_),
          orbital_f_dn(orbital_f_dn_),
          d1_exponent(d1_exponent_),
          d2_exponent(d2_exponent_),
          f_up_exponent(f_up_exponent_),
          f_dn_exponent(f_dn_exponent_)
        {
        }

    class Amplitude : public Wavefunction<AmplitudeType>::Amplitude
    {
    private:
        CeperleyMatrix<AmplitudeType> m_cmat_d1, m_cmat_d2, m_cmat_f_up, m_cmat_f_dn;
        int m_partial_update_step;

        // the following variables only need be set when a move is in progress
        unsigned int m_up_particles_in_progress, m_down_particles_in_progress;
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const DMetalWavefunction> &wf, const PositionArguments &r_);

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

    virtual boost::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> create_nonzero_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction<AmplitudeType> > &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const override;

    virtual boost::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction<AmplitudeType> > &this_ptr, const PositionArguments &r) const override
        {
            BOOST_ASSERT(this == this_ptr.get());
            return boost::make_shared<Amplitude>(boost::dynamic_pointer_cast<const DMetalWavefunction>(this_ptr), r);
        }

    virtual unsigned int get_N_species (void) const override
        {
            return 2;
        }

    virtual unsigned int get_N_filled (unsigned int species) const override
        {
            BOOST_ASSERT(species < 2);
            return (species == 0 ? orbital_d1 : orbital_d2)->get_N_filled();
        }
};

#define VMC_SUPPORTED_TYPE(type) extern template class DMetalWavefunction<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_TYPE

#endif
