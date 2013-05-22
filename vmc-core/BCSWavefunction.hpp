#ifndef _VMC_BCS_WAVEFUNCTION_HPP
#define _VMC_BCS_WAVEFUNCTION_HPP

#include <vector>

#include <boost/make_shared.hpp>
#include <boost/assert.hpp>

// we always want to include vmc-typedefs.hpp before including Eigen
#include "vmc-typedefs.hpp"
#include <Eigen/Dense>

#include "Wavefunction.hpp"
#include "CeperleyMatrix.hpp"
#include "JastrowFactor.hpp"

/**
 * Projected BCS wave function
 *
 * Assumes unpolarized state.
 */
template <typename AmplitudeType>
class BCSWavefunction : public Wavefunction
{
public:
    const Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> phi;
    const unsigned int M; // number of particles of each species
    const boost::shared_ptr<const JastrowFactor> jastrow;

    BCSWavefunction (const boost::shared_ptr<const Lattice> &lattice_, const Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> &phi_, unsigned int M_, const boost::shared_ptr<const JastrowFactor> &jastrow_=boost::shared_ptr<const JastrowFactor>())
        : Wavefunction(lattice_),
          phi(phi_),
          M(M_),
          jastrow(jastrow_)
        {
            BOOST_ASSERT(M > 0);
            BOOST_ASSERT(2 * M <= lattice->total_sites());
            BOOST_ASSERT(lattice->total_sites() == phi.rows());
            BOOST_ASSERT(lattice->total_sites() == phi.cols());
        }

    class Amplitude : public Wavefunction::Amplitude
    {
    private:
        CeperleyMatrix<AmplitudeType> m_cmat;
        Big<AmplitudeType> m_current_jastrow, m_old_jastrow;
        unsigned int m_partial_update_step;
        Move m_current_move;

    public:
        Amplitude (const boost::shared_ptr<const BCSWavefunction> &wf_, const PositionArguments &r_);

        virtual Move propose_move (RandomNumberGenerator &rng) const;

    private:
        virtual void perform_move_ (const Move &move) override;

        template <bool first_pass>
        void do_perform_move (const Move &move);

        virtual Big<AmplitudeType> psi_ (void) const override;

        virtual void finish_move_ (void) override;

        virtual void cancel_move_ (void) override;

        virtual void swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species) override;

        virtual void reset_ (const PositionArguments &r_) override;

        virtual boost::shared_ptr<Wavefunction::Amplitude> clone_ (void) const override;

        void reinitialize (void);

        virtual void check_for_numerical_error (void) const override;
    };

    virtual boost::shared_ptr<Wavefunction::Amplitude> create_nonzero_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const override;

    virtual boost::shared_ptr<Wavefunction::Amplitude> create_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, const PositionArguments &r) const override
        {
            BOOST_ASSERT(this == this_ptr.get());
            return boost::make_shared<Amplitude>(boost::dynamic_pointer_cast<const BCSWavefunction>(this_ptr), r);
        }

    virtual unsigned int get_N_species (void) const override
        {
            return 2;
        }

    virtual unsigned int get_N_filled (unsigned int species) const override
        {
            BOOST_ASSERT(species < 2);
            // assumes half filling
            return M;
        }

    // CYTHON-LIMITATION: http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#c-left-values
    static inline void set_matrix_coeff (Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> &mat, unsigned int row, unsigned int col, AmplitudeType value)
        {
            mat(row, col) = value;
        }
};

#define VMC_SUPPORTED_TYPE(type) extern template class BCSWavefunction<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_TYPE

#endif
