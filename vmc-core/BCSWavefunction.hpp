#ifndef _BCS_WAVEFUNCTION_HPP
#define _BCS_WAVEFUNCTION_HPP

#include <vector>

#include <boost/make_shared.hpp>

#include "Wavefunction.hpp"
#include "CeperleyMatrix.hpp"

/**
 * Projected BCS wave function at half-filling
 */
class BCSWavefunction : public Wavefunction
{
public:
    const std::vector<complex_t> phi;

    BCSWavefunction (const boost::shared_ptr<const Lattice> &lattice_, const std::vector<complex_t> &phi_)
        : Wavefunction(lattice_),
          phi(phi_)
        {
        }

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        // Because we don't yet have a way of doing row-column updates in a
        // single step, on each step we copy m_cmat to m_new_cmat, update it
        // twice, and copy it back if we choose to accept the move.
        bool m_update_in_progress;
        CeperleyMatrix<amplitude_t> m_cmat;

    public:
        Amplitude (const boost::shared_ptr<const BCSWavefunction> &wf_, const PositionArguments &r_);

        virtual Move propose_move (RandomNumberGenerator &rng) const;

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
            return lattice->total_sites() / 2;
        }
};

#endif
