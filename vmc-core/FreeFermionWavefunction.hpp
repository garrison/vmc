#ifndef _VMC_FREE_FERMION_WAVEFUNCTION_HPP
#define _VMC_FREE_FERMION_WAVEFUNCTION_HPP

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

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
class FreeFermionWavefunction : public Wavefunction
{
public:
    const std::vector<boost::shared_ptr<const OrbitalDefinitions> > orbital_def;
    const boost::shared_ptr<const JastrowFactor> jastrow;

    FreeFermionWavefunction (const std::vector<boost::shared_ptr<const OrbitalDefinitions> > &orbital_def_, const boost::shared_ptr<const JastrowFactor> &jastrow_=boost::shared_ptr<const JastrowFactor>());

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        std::vector<CeperleyMatrix<amplitude_t> > m_cmat;
        Big<amplitude_t> m_current_jastrow, m_old_jastrow;
        unsigned int m_partial_update_step;
        std::vector<bool> m_species_move_in_progress;
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_);

    private:
        virtual void perform_move_ (const Move &move) override;

        template <bool first_pass>
        void do_perform_move (const Move &move);

        virtual Big<amplitude_t> psi_ (void) const override;

        virtual void finish_move_ (void) override;

        virtual void cancel_move_ (void) override;

        virtual void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species) override;

        virtual void reset_ (const PositionArguments &r_) override;

        virtual boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const override;

        void reinitialize (void);

        virtual void check_for_numerical_error (void) const override;

        unsigned int get_N_species (void) const
            {
                return m_species_move_in_progress.size();
            }
    };

    virtual boost::shared_ptr<Wavefunction::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, const PositionArguments &r) const override
        {
            BOOST_ASSERT(this == this_ptr.get());
            return boost::make_shared<Amplitude>(boost::dynamic_pointer_cast<const FreeFermionWavefunction>(this_ptr), r);
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

#endif
