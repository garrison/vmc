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
        real_t m_current_jastrow, m_old_jastrow;
        unsigned int m_partial_update_step;
        std::vector<bool> m_species_move_in_progress;
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_);

    private:
        void perform_move_ (const Move &move);

        template <bool first_pass>
        void do_perform_move (const Move &move);

        Big<amplitude_t> psi_ (void) const;

        void finish_move_ (void);

        void cancel_move_ (void);

        void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species);

        void reset_ (const PositionArguments &r_);

        boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const;

        void reinitialize (void);

        unsigned int get_N_species (void) const
            {
                return m_species_move_in_progress.size();
            }
    };

    boost::shared_ptr<Wavefunction::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, const PositionArguments &r) const
        {
            BOOST_ASSERT(this == this_ptr.get());
            return boost::make_shared<Amplitude>(boost::dynamic_pointer_cast<const FreeFermionWavefunction>(this_ptr), r);
        }

    unsigned int get_N_species (void) const
        {
            return orbital_def.size();
        }

    unsigned int get_N_filled (unsigned int species) const
        {
            BOOST_ASSERT(species < get_N_species());
            return orbital_def[species]->get_N_filled();
        }
};

#endif
