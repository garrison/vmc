#include <vector>

#include <boost/assert.hpp>
#include <boost/cast.hpp>

#include "DMetalWavefunction.hpp"
#include "vmc-math-utils.hpp"
#include "random-configuration.hpp"

template <typename AmplitudeType>
DMetalWavefunction<AmplitudeType>::Amplitude::Amplitude (const std::shared_ptr<const DMetalWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction<AmplitudeType>::Amplitude(wf_, r_),
      m_partial_update_step(0)
{
    reinitialize();
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::perform_move_ (const Move &move)
{
    // we require that m_partial_update_step == 0 between moves; otherwise, psi_() will
    // return zero when it shouldn't.
    BOOST_ASSERT(m_partial_update_step == 0);

    m_down_particles_in_progress = 0;
    for (unsigned int i = 0; i < move.size(); ++i) {
        BOOST_ASSERT(move[i].particle.species < 2);
        m_down_particles_in_progress += move[i].particle.species;
    }
    m_up_particles_in_progress = move.size() - m_down_particles_in_progress;

    m_current_move = move;

    do_perform_move<true>(move);
}

// During perform_move(), we want to short-circuit out as soon as any of the
// determinants is zero, as we already know what psi_() will be.  However,
// doing so means we must later make a second pass to finish all the
// determinant updates if finish_update() is called.  This templated function
// allows us to use the same code in both passes, with only minor differences.
template <typename AmplitudeType>
template <bool first_pass>
void DMetalWavefunction<AmplitudeType>::Amplitude::do_perform_move (const Move &move)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(this->wf.get());

    const unsigned int M = wf_->orbital_f_up->get_N_filled();

    if (first_pass) {
        // explicitly enforce Gutzwiller projection
        for (unsigned int i = 0; i < move.size(); ++i) {
            if (this->r.is_occupied(move[i].destination, move[i].particle.species ^ 1)) {
                m_partial_update_step = 4;
                return;
            }
        }
    }

    switch (first_pass ? 4 : m_partial_update_step) {
    case 4:
        {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> d2_cols;
            for (unsigned int i = 0; i < move.size(); ++i) {
                const Particle &particle = move[i].particle;
                const unsigned int col_index = (particle.species == 0) ? particle.index : particle.index + M;
                d2_cols.push_back(std::make_pair(col_index, move[i].destination));
            }
            m_cmat_d2.update_columns(d2_cols, wf_->orbital_d2->get_orbitals());
        }
        if (first_pass && m_cmat_d2.is_singular()) {
            m_partial_update_step = 3;
            return;
        }

    case 3:
        {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> d1_cols;
            for (unsigned int i = 0; i < move.size(); ++i) {
                const Particle &particle = move[i].particle;
                const unsigned int col_index = (particle.species == 0) ? particle.index : particle.index + M;
                d1_cols.push_back(std::make_pair(col_index, move[i].destination));
            }
            m_cmat_d1.update_columns(d1_cols, wf_->orbital_d1->get_orbitals());
        }
        if (first_pass && m_cmat_d1.is_singular()) {
            m_partial_update_step = 2;
            return;
        }

    case 2:
        if (m_up_particles_in_progress) {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> f_up_cols;
            for (unsigned int i = 0; i < move.size(); ++i) {
                if (move[i].particle.species == 0)
                    f_up_cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
            }
            BOOST_ASSERT(f_up_cols.size() == m_up_particles_in_progress);
            if (m_up_particles_in_progress)
                m_cmat_f_up.update_columns(f_up_cols, wf_->orbital_f_up->get_orbitals());
        }
        if (first_pass && m_cmat_f_up.is_singular()) {
            m_partial_update_step = 1;
            return;
        }

    case 1:
        if (m_down_particles_in_progress) {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> f_dn_cols;
            for (unsigned int i = 0; i < move.size(); ++i) {
                if (move[i].particle.species != 0)
                    f_dn_cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
            }
            BOOST_ASSERT(f_dn_cols.size() == m_down_particles_in_progress);
            if (m_down_particles_in_progress)
                m_cmat_f_dn.update_columns(f_dn_cols, wf_->orbital_f_dn->get_orbitals());
        }
        m_partial_update_step = 0;

    case 0: ;
    }
}

template <typename AmplitudeType>
Big<AmplitudeType> DMetalWavefunction<AmplitudeType>::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return Big<AmplitudeType>();

    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(this->wf.get());

    return (complex_pow(m_cmat_d1.get_determinant(), wf_->d1_exponent)
            * complex_pow(m_cmat_d2.get_determinant(), wf_->d2_exponent)
            * complex_pow(m_cmat_f_up.get_determinant(), wf_->f_up_exponent)
            * complex_pow(m_cmat_f_dn.get_determinant(), wf_->f_dn_exponent));
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    m_cmat_d1.finish_columns_update();
    m_cmat_d2.finish_columns_update();
    if (m_up_particles_in_progress)
        m_cmat_f_up.finish_columns_update();
    if (m_down_particles_in_progress)
        m_cmat_f_dn.finish_columns_update();
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::cancel_move_ (void)
{
    switch (m_partial_update_step) {
    case 0:
        if (m_down_particles_in_progress)
            m_cmat_f_dn.cancel_columns_update();
    case 1:
        if (m_up_particles_in_progress)
            m_cmat_f_up.cancel_columns_update();
    case 2:
        m_cmat_d1.cancel_columns_update();
    case 3:
        m_cmat_d2.cancel_columns_update();
    case 4:
        ;
    }

    m_partial_update_step = 0;
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(this->wf.get());
    const unsigned int M = wf_->orbital_f_up->get_N_filled();
    const unsigned int particle1_column_index = (species == 0) ? particle1_index : particle1_index + M;
    const unsigned int particle2_column_index = (species == 0) ? particle2_index : particle2_index + M;
    m_cmat_d1.swap_columns(particle1_column_index, particle2_column_index);
    m_cmat_d2.swap_columns(particle1_column_index, particle2_column_index);
    (species == 0 ? m_cmat_f_up : m_cmat_f_dn).swap_columns(particle1_index, particle2_index);
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::reset_ (const PositionArguments &r_)
{
    this->r = r_;
    reinitialize();
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::reinitialize (void)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(this->wf.get());

    BOOST_ASSERT(this->r.get_N_species() == 2);

    BOOST_ASSERT(this->r.get_N_sites() == wf_->orbital_d1->get_N_sites());
    BOOST_ASSERT(wf_->orbital_d1->get_lattice_ptr() == wf_->orbital_d2->get_lattice_ptr());
    BOOST_ASSERT(wf_->orbital_d1->get_lattice_ptr() == wf_->orbital_f_up->get_lattice_ptr());
    BOOST_ASSERT(wf_->orbital_d1->get_lattice_ptr() == wf_->orbital_f_dn->get_lattice_ptr());

    BOOST_ASSERT(this->r.get_N_filled(0) + this->r.get_N_filled(1) == wf_->orbital_d1->get_N_filled());
    BOOST_ASSERT(this->r.get_N_filled(0) + this->r.get_N_filled(1) == wf_->orbital_d2->get_N_filled());
    BOOST_ASSERT(this->r.get_N_filled(0) == wf_->orbital_f_up->get_N_filled());
    BOOST_ASSERT(this->r.get_N_filled(1) == wf_->orbital_f_dn->get_N_filled());

    const unsigned int N = wf_->orbital_d1->get_N_filled();
    const unsigned int M = wf_->orbital_f_up->get_N_filled();

    Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> mat_d1(N, N);
    Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> mat_d2(N, N);
    Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> mat_f_up(M, M);
    Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> mat_f_dn(N - M, N - M);

    for (unsigned int i = 0; i < this->r.get_N_filled(0); ++i) {
        const Particle particle(i, 0);
        mat_d1.col(i) = wf_->orbital_d1->at_position(this->r[particle]);
        mat_d2.col(i) = wf_->orbital_d2->at_position(this->r[particle]);
        mat_f_up.col(i) = wf_->orbital_f_up->at_position(this->r[particle]);
    }

    for (unsigned int i = 0; i < this->r.get_N_filled(1); ++i) {
        const Particle particle(i, 1);
        mat_d1.col(i + M) = wf_->orbital_d1->at_position(this->r[particle]);
        mat_d2.col(i + M) = wf_->orbital_d2->at_position(this->r[particle]);
        mat_f_dn.col(i) = wf_->orbital_f_dn->at_position(this->r[particle]);
    }

    m_cmat_d1 = CeperleyMatrix<AmplitudeType>(mat_d1);
    m_cmat_d2 = CeperleyMatrix<AmplitudeType>(mat_d2);
    m_cmat_f_up = CeperleyMatrix<AmplitudeType>(mat_f_up);
    m_cmat_f_dn = CeperleyMatrix<AmplitudeType>(mat_f_dn);
}

template <typename AmplitudeType>
void DMetalWavefunction<AmplitudeType>::Amplitude::check_for_numerical_error (void) const
{
    m_cmat_d1.check_for_numerical_error();
    m_cmat_d2.check_for_numerical_error();
    m_cmat_f_up.check_for_numerical_error();
    m_cmat_f_dn.check_for_numerical_error();
}

template <typename AmplitudeType>
std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> DMetalWavefunction<AmplitudeType>::Amplitude::clone_ (void) const
{
    return std::make_shared<DMetalWavefunction<AmplitudeType>::Amplitude>(*this);
}

template <typename AmplitudeType>
std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> DMetalWavefunction<AmplitudeType>::create_nonzero_wavefunctionamplitude (const std::shared_ptr<const Wavefunction<AmplitudeType> > &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const
{
    const unsigned int N = orbital_d1->get_N_filled();
    const unsigned int M = orbital_f_up->get_N_filled();

    // take into account the Gutzwiller projection
    while (n_attempts--) {
        std::vector<std::vector<unsigned int> > vv(2);
        vv[0] = some_random_configuration(N, *this->lattice, rng);
        for (unsigned int i = M; i < N; ++i)
            vv[1].push_back(vv[0][i]);
        vv[0].resize(M);
        std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> wfa(create_wavefunctionamplitude(this_ptr, PositionArguments(vv, this->lattice->total_sites())));
        if (wfa->is_nonzero())
            return wfa;
    }
    return std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude>();
}

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) template class DMetalWavefunction<type>
#include "vmc-supported-types.hpp"
