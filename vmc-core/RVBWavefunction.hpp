#ifndef _RVB_WAVEFUNCTION_HPP
#define _RVB_WAVEFUNCTION_HPP

#include <vector>

#include "Wavefunction.hpp"
#include "CeperleyMatrix.hpp"

/**
 * RVB (projected BCS wave function at half-filling)
 */
class RVBWavefunction : public Wavefunction
{
public:
    const std::vector<complex_t> phi;

    RVBWavefunction (const boost::shared_ptr<const Lattice> &lattice_, const std::vector<complex_t> &phi_)
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
        CeperleyMatrix<amplitude_t> m_cmat, m_new_cmat;

    public:
        Amplitude (const boost::shared_ptr<const RVBWavefunction> &wf_, const PositionArguments &r_);

        virtual Move propose_move (RandomNumberGenerator &rng) const;

    private:
        void perform_move_ (const Move &move);

        amplitude_t psi_ (void) const;

        void finish_move_ (void);

        void cancel_move_ (void);

        void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species);

        void reset_ (const PositionArguments &r_);

        boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const;

        void reset_with_random_configuration (RandomNumberGenerator &rng);

        void reinitialize (void);
    };
};

#endif
