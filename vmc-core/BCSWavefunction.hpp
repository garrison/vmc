#ifndef _BCS_WAVEFUNCTION_HPP
#define _BCS_WAVEFUNCTION_HPP

#include <vector>

#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

#include "Wavefunction.hpp"
#include "CeperleyMatrix.hpp"
#include "JastrowFactor.hpp"

/**
 * Projected BCS wave function
 *
 * Assumes unpolarized state.
 *
 * Assumes periodic boundary conditions.
 */
class BCSWavefunction : public Wavefunction
{
public:
    const std::vector<complex_t> phi;
    const unsigned int M; // number of particles of each species
    const boost::shared_ptr<const JastrowFactor> jastrow;

    BCSWavefunction (const boost::shared_ptr<const Lattice> &lattice_, const std::vector<complex_t> &phi_, unsigned int M_, const boost::shared_ptr<const JastrowFactor> &jastrow_=boost::shared_ptr<const JastrowFactor>())
        : Wavefunction(lattice_),
          phi(phi_),
          M(M_),
          jastrow(jastrow_)
        {
            BOOST_ASSERT(M > 0);
            BOOST_ASSERT(2 * M <= lattice->total_sites());
            BOOST_ASSERT(lattice->total_sites() == phi.size());
        }

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        CeperleyMatrix<amplitude_t> m_cmat;
        real_t m_current_jastrow, m_old_jastrow;
        unsigned int m_partial_update_step;
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const BCSWavefunction> &wf_, const PositionArguments &r_);

        virtual Move propose_move (RandomNumberGenerator &rng) const;

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

    boost::shared_ptr<Wavefunction::Amplitude> create_nonzero_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const;

    boost::shared_ptr<Wavefunction::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, const PositionArguments &r) const
        {
            BOOST_ASSERT(this == this_ptr.get());
            return boost::make_shared<Amplitude>(boost::shared_polymorphic_downcast<const BCSWavefunction>(this_ptr), r);
        }

    unsigned int get_N_species (void) const
        {
            return 2;
        }

    unsigned int get_N_filled (unsigned int species) const
        {
            BOOST_ASSERT(species < 2);
            // assumes half filling
            return M;
        }
};

#endif
