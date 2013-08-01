#ifndef _VMC_FREE_FERMION_WAVEFUNCTION_HPP
#define _VMC_FREE_FERMION_WAVEFUNCTION_HPP

#include <vector>

#include "Wavefunction.hpp"
#include "PositionArguments.hpp"
#include "CeperleyMatrix.hpp"
#include "OrbitalDefinitions.hpp"
#include "JastrowFactor.hpp"

/**
 * Free fermion wave function
 *
 * (a single determinant for each species, with no Jastrow factor)
 */
template <typename AmplitudeType>
class FreeFermionWavefunction : public Wavefunction<AmplitudeType>
{
public:
    const std::vector<std::shared_ptr<const OrbitalDefinitions<AmplitudeType> > > orbital_def;
    const std::shared_ptr<const JastrowFactor<AmplitudeType> > jastrow;

    FreeFermionWavefunction (const std::vector<std::shared_ptr<const OrbitalDefinitions<AmplitudeType> > > &orbital_def_, const std::shared_ptr<const JastrowFactor<AmplitudeType> > &jastrow_=std::shared_ptr<const JastrowFactor<AmplitudeType> >());

    class Amplitude : public Wavefunction<AmplitudeType>::Amplitude
    {
    private:
        std::vector<CeperleyMatrix<AmplitudeType> > m_cmat;
        Big<AmplitudeType> m_current_jastrow, m_old_jastrow;
        unsigned int m_partial_update_step;
        std::vector<bool> m_species_move_in_progress;
        Move m_current_move;

    public:
        Amplitude (const std::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_);

    private:
        virtual void perform_move_ (const Move &move) override;

        template <bool first_pass>
        void do_perform_move (const Move &move);

        virtual Big<AmplitudeType> psi_ (void) const override;

        virtual void finish_move_ (void) override;

        virtual void cancel_move_ (void) override;

        virtual void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species) override;

        virtual void reset_ (const PositionArguments &r_) override;

        virtual std::unique_ptr<typename Wavefunction<AmplitudeType>::Amplitude> clone_ (void) const override;

        void reinitialize (void);

        virtual void check_for_numerical_error (void) const override;

        unsigned int get_N_species (void) const
            {
                return m_species_move_in_progress.size();
            }
    };

    virtual std::unique_ptr<typename Wavefunction<AmplitudeType>::Amplitude> create_wavefunctionamplitude (const std::shared_ptr<const Wavefunction<AmplitudeType> > &this_ptr, const PositionArguments &r) const override
        {
            BOOST_ASSERT(this == this_ptr.get());
            return std::make_unique<Amplitude>(std::dynamic_pointer_cast<const FreeFermionWavefunction>(this_ptr), r);
        }

    virtual unsigned int get_N_species (void) const override
        {
            return orbital_def.size();
        }

    virtual unsigned int get_N_filled (unsigned int species) const override
        {
            BOOST_ASSERT(species < get_N_species());
            return orbital_def[species]->get_N_filled();
        }
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class FreeFermionWavefunction<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif
