#include <vector>

#include <boost/assert.hpp>
#include <boost/make_shared.hpp>
#include <boost/cast.hpp>

#include "vmc-math-utils.hpp"
#include "DMetalWavefunction.hpp"
#include "random-configuration.hpp"

DMetalWavefunction::Amplitude::Amplitude (const boost::shared_ptr<const DMetalWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction::Amplitude(wf_, r_),
      m_partial_update_step(0)
{
    reinitialize();
}

void DMetalWavefunction::Amplitude::perform_move_ (const Move &move)
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
template <bool first_pass>
void DMetalWavefunction::Amplitude::do_perform_move (const Move &move)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(wf.get());

    const unsigned int M = wf_->orbital_f_up->get_N_filled();

    if (first_pass) {
        // explicitly enforce Gutzwiller projection
        for (unsigned int i = 0; i < move.size(); ++i) {
            if (r.is_occupied(move[i].destination, move[i].particle.species ^ 1)) {
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
        if (first_pass && m_cmat_d2.get_determinant() == amplitude_t(0)) {
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
        if (first_pass && m_cmat_d1.get_determinant() == amplitude_t(0)) {
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
        if (first_pass && m_cmat_f_up.get_determinant() == amplitude_t(0)) {
            m_partial_update_step = 1;
            return;
        }

    case 1:
        if (m_down_particles_in_progress) {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> f_down_cols;
            for (unsigned int i = 0; i < move.size(); ++i) {
                if (move[i].particle.species != 0)
                    f_down_cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
            }
            BOOST_ASSERT(f_down_cols.size() == m_down_particles_in_progress);
            if (m_down_particles_in_progress)
                m_cmat_f_down.update_columns(f_down_cols, wf_->orbital_f_down->get_orbitals());
        }
        m_partial_update_step = 0;

    case 0: ;
    }
}

amplitude_t DMetalWavefunction::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return amplitude_t(0);

    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(wf.get());

    return (complex_pow(m_cmat_d1.get_determinant(), wf_->d1_exponent)
            * complex_pow(m_cmat_d2.get_determinant(), wf_->d2_exponent)
            * complex_pow(m_cmat_f_up.get_determinant(), wf_->f_up_exponent)
            * complex_pow(m_cmat_f_down.get_determinant(), wf_->f_down_exponent));
}

void DMetalWavefunction::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    m_cmat_d1.finish_columns_update();
    m_cmat_d2.finish_columns_update();
    if (m_up_particles_in_progress)
        m_cmat_f_up.finish_columns_update();
    if (m_down_particles_in_progress)
        m_cmat_f_down.finish_columns_update();
}

void DMetalWavefunction::Amplitude::cancel_move_ (void)
{
    switch (m_partial_update_step) {
    case 0:
        if (m_down_particles_in_progress)
            m_cmat_f_down.cancel_columns_update();
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

void DMetalWavefunction::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(wf.get());
    const unsigned int M = wf_->orbital_f_up->get_N_filled();
    const unsigned int particle1_column_index = (species == 0) ? particle1_index : particle1_index + M;
    const unsigned int particle2_column_index = (species == 0) ? particle2_index : particle2_index + M;
    m_cmat_d1.swap_columns(particle1_column_index, particle2_column_index);
    m_cmat_d2.swap_columns(particle1_column_index, particle2_column_index);
    (species == 0 ? m_cmat_f_up : m_cmat_f_down).swap_columns(particle1_index, particle2_index);
}

void DMetalWavefunction::Amplitude::reset_ (const PositionArguments &r_)
{
    r = r_;
    reinitialize();
}

void DMetalWavefunction::Amplitude::reinitialize (void)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(wf.get());

    BOOST_ASSERT(r.get_N_species() == 2);

    BOOST_ASSERT(r.get_N_sites() == wf_->orbital_d1->get_N_sites());
    BOOST_ASSERT(wf_->orbital_d1->get_lattice_ptr() == wf_->orbital_d2->get_lattice_ptr());
    BOOST_ASSERT(wf_->orbital_d1->get_lattice_ptr() == wf_->orbital_f_up->get_lattice_ptr());
    BOOST_ASSERT(wf_->orbital_d1->get_lattice_ptr() == wf_->orbital_f_down->get_lattice_ptr());

    BOOST_ASSERT(r.get_N_filled(0) + r.get_N_filled(1) == wf_->orbital_d1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) + r.get_N_filled(1) == wf_->orbital_d2->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) == wf_->orbital_f_up->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(1) == wf_->orbital_f_down->get_N_filled());

    const unsigned int N = wf_->orbital_d1->get_N_filled();
    const unsigned int M = wf_->orbital_f_up->get_N_filled();

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_d1(N, N);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_d2(N, N);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_f_up(M, M);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_f_down(N - M, N - M);

    for (unsigned int i = 0; i < r.get_N_filled(0); ++i) {
        const Particle particle(i, 0);
        mat_d1.col(i) = wf_->orbital_d1->at_position(r[particle]);
        mat_d2.col(i) = wf_->orbital_d2->at_position(r[particle]);
        mat_f_up.col(i) = wf_->orbital_f_up->at_position(r[particle]);
    }

    for (unsigned int i = 0; i < r.get_N_filled(1); ++i) {
        const Particle particle(i, 1);
        mat_d1.col(i + M) = wf_->orbital_d1->at_position(r[particle]);
        mat_d2.col(i + M) = wf_->orbital_d2->at_position(r[particle]);
        mat_f_down.col(i) = wf_->orbital_f_down->at_position(r[particle]);
    }

    m_cmat_d1 = mat_d1;
    m_cmat_d2 = mat_d2;
    m_cmat_f_up = mat_f_up;
    m_cmat_f_down = mat_f_down;
}

boost::shared_ptr<Wavefunction::Amplitude> DMetalWavefunction::Amplitude::clone_ (void) const
{
    return boost::make_shared<DMetalWavefunction::Amplitude>(*this);
}

void DMetalWavefunction::Amplitude::reset_with_random_configuration (RandomNumberGenerator &rng)
{
    const DMetalWavefunction *wf_ = boost::polymorphic_downcast<const DMetalWavefunction *>(wf.get());

    const unsigned int N = wf_->orbital_d1->get_N_filled();
    const unsigned int M = wf_->orbital_f_up->get_N_filled();

    // take into account the Gutzwiller projection
    std::vector<std::vector<unsigned int> > vv(2);
    vv[0] = some_random_configuration(N, *wf->lattice, rng);
    for (unsigned int i = M; i < N; ++i)
        vv[1].push_back(vv[0][i]);
    vv[0].resize(M);

    reset(PositionArguments(vv, wf->lattice->total_sites()));
}
