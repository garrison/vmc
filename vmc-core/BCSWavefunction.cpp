#include <vector>
#include <algorithm>

#include <boost/cast.hpp>

#include "BCSWavefunction.hpp"
#include "RandomNumberGenerator.hpp"
#include "random-combination.hpp"
#include "random-configuration.hpp"
#include "random-move.hpp"

template <typename AmplitudeType>
BCSWavefunction<AmplitudeType>::Amplitude::Amplitude (const std::shared_ptr<const BCSWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction<AmplitudeType>::Amplitude(wf_, r_),
      m_current_jastrow(1),
      m_partial_update_step(0)
{
    reinitialize();
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::perform_move_ (const Move &move)
{
    // we require that m_partial_update_step == 0 between moves; otherwise,
    // psi_() will return zero when it shouldn't.
    BOOST_ASSERT(m_partial_update_step == 0);

    m_current_move = move;

    m_old_jastrow = m_current_jastrow;

    do_perform_move<true>(move);
}

// During perform_move(), we want to short-circuit out as soon as any of the
// determinants is zero, as we already know what psi_() will be.  However,
// doing so means we must later make a second pass to finish all the
// determinant updates if finish_update() is called.  This templated function
// allows us to use the same code in both passes, with only minor differences.
template <typename AmplitudeType>
template <bool first_pass>
void BCSWavefunction<AmplitudeType>::Amplitude::do_perform_move (const Move &move)
{
    const BCSWavefunction *wf_ = boost::polymorphic_downcast<const BCSWavefunction *>(this->wf.get());
    const unsigned int M = wf_->M;

    if (first_pass) {
        // explicitly enforce Gutzwiller projection
        for (unsigned int i = 0; i < move.size(); ++i) {
            if (this->r.is_occupied(move[i].destination, move[i].particle.species ^ 1)) {
                m_partial_update_step = 2;
                return;
            }
        }
    }

    switch (first_pass ? 2 : m_partial_update_step) {
    case 2:
        if (wf_->jastrow) {
            m_current_jastrow = wf_->jastrow->compute_jastrow(this->r);
            if (m_current_jastrow.is_zero()) {
                m_partial_update_step = 1;
                return;
            }
        } else {
            BOOST_ASSERT(m_current_jastrow.get_base() == 1. && m_current_jastrow.get_exponent() == 0.);
        }

    case 1:
        {
            Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> srcmat(M, M);
            lw_vector<unsigned int, MAX_MOVE_SIZE> rows, cols;

            for (unsigned int j = 0; j < move.size(); ++j) {
                const unsigned int particle_index = move[j].particle.index;
                const unsigned int destination = move[j].destination;

                if (move[j].particle.species == 0) {
                    const std::vector<unsigned int> & dn_pos = this->r.r_vector(1);
                    for (unsigned int i = 0; i < M; ++i)
                        srcmat(particle_index, i) = wf_->phi(destination, dn_pos[i]);
                    rows.push_back(particle_index);
                } else {
                    BOOST_ASSERT(move[j].particle.species == 1);
                    const std::vector<unsigned int> & up_pos = this->r.r_vector(0);
                    for (unsigned int i = 0; i < M; ++i)
                        srcmat(i, particle_index) = wf_->phi(up_pos[i], destination);
                    cols.push_back(particle_index);
                }
            }

            m_cmat.update_rows_and_columns(rows, cols, srcmat);
        }

        m_partial_update_step = 0;

    case 0: ;
    }
}

template <typename AmplitudeType>
Big<AmplitudeType> BCSWavefunction<AmplitudeType>::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return Big<AmplitudeType>();

    return m_cmat.get_determinant() * m_current_jastrow;
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    m_cmat.finish_rows_and_columns_update();
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::cancel_move_ (void)
{
    if (m_partial_update_step == 0)
        m_cmat.cancel_rows_and_columns_update();
    m_current_jastrow = m_old_jastrow; // this is always safe as a way to reset
                                       // m_current_jastrow, even during a
                                       // partial update
    m_partial_update_step = 0;
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    if (species == 0) {
        m_cmat.swap_rows(particle1_index, particle2_index);
    } else {
        BOOST_ASSERT(species == 1);
        m_cmat.swap_columns(particle1_index, particle2_index);
    }
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::reset_ (const PositionArguments &r_)
{
    this->r = r_;
    reinitialize();
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::reinitialize (void)
{
    const BCSWavefunction *wf_ = boost::polymorphic_downcast<const BCSWavefunction *>(this->wf.get());
    const std::shared_ptr<const Lattice> &lattice = this->wf->lattice;
    const unsigned int M = wf_->M;

    BOOST_ASSERT(this->r.get_N_sites() == lattice->total_sites());
    BOOST_ASSERT(this->r.get_N_species() == 2);
    BOOST_ASSERT(this->r.get_N_filled(0) == M);
    BOOST_ASSERT(this->r.get_N_filled(1) == M);

    if (wf_->jastrow) {
        m_current_jastrow = wf_->jastrow->compute_jastrow(this->r);
    }

    Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> mat_phi(M, M);
    const std::vector<unsigned int> & up_pos = this->r.r_vector(0);
    const std::vector<unsigned int> & dn_pos = this->r.r_vector(1);
    for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < M; ++j) {
            mat_phi(i, j) = wf_->phi(up_pos[i], dn_pos[j]);
        }
    }

    m_cmat = CeperleyMatrix<AmplitudeType>(mat_phi);
}

template <typename AmplitudeType>
void BCSWavefunction<AmplitudeType>::Amplitude::check_for_numerical_error (void) const
{
    m_cmat.check_for_numerical_error();
}

template <typename AmplitudeType>
std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> BCSWavefunction<AmplitudeType>::Amplitude::clone_ (void) const
{
    return std::make_shared<BCSWavefunction<AmplitudeType>::Amplitude>(*this);
}

template <typename AmplitudeType>
std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> BCSWavefunction<AmplitudeType>::create_nonzero_wavefunctionamplitude (const std::shared_ptr<const Wavefunction<AmplitudeType> > &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const
{
    const unsigned int M = get_N_filled(0);
    const unsigned int N = this->lattice->total_sites();

    BOOST_ASSERT(2 * M <= N);

    while (n_attempts--) {
        // Take into account Gutzwiller projection [at most one spinon ("particle") per site]
        std::vector<std::vector<unsigned int> > vv(2);
        vv[0] = some_random_configuration(M, *this->lattice, rng);
        std::vector<unsigned int> occupied_sites(N);
        for (unsigned int i = 0; i < vv[0].size(); ++i) {
            BOOST_ASSERT(vv[0][i] < N);
            BOOST_ASSERT(occupied_sites[vv[0][i]] == 0);
            ++occupied_sites[vv[0][i]];
        }
        std::vector<unsigned int> unoccupied_site_indices; // unoccupied by the first species, that is
        unoccupied_site_indices.reserve(N - M);
        for (unsigned int i = 0; i < N; ++i) {
            if (occupied_sites[i] == 0)
                unoccupied_site_indices.push_back(i);
        }
        BOOST_ASSERT(unoccupied_site_indices.size() == N - M);
        std::vector<unsigned int> choices;
        random_combination(choices, M, N - M, rng);
        for (unsigned int i = 0; i < choices.size(); ++i)
            vv[1].push_back(unoccupied_site_indices[i]);

        BOOST_ASSERT(vv[0].size() == M);
        BOOST_ASSERT(vv[1].size() == M);

        std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> wfa(create_wavefunctionamplitude(this_ptr, PositionArguments(vv, this->lattice->total_sites())));
        if (wfa->is_nonzero())
            return wfa;
    }
    return std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude>();
}

template <typename AmplitudeType>
Move BCSWavefunction<AmplitudeType>::Amplitude::propose_move (RandomNumberGenerator &rng) const
{
    Move move;

    // choose a particle and plan to move it to a site that doesn't already
    // contain a particle of the same spin
    const Particle particle(choose_random_particle(this->r, rng));
    const unsigned int target_site_index = plan_particle_move_to_nearby_empty_site(particle, this->r, *this->wf->lattice, rng);
    if (target_site_index == this->r[particle])
        return Move();
    move.push_back(SingleParticleMove(particle, target_site_index));

    // if the target site has a particle of opposite spin, we want it to swap
    // positions with our original particle
    BOOST_ASSERT(this->wf->get_N_species() == 2);
    const unsigned int other_species = particle.species ^ 1;
    if (this->r.is_occupied(target_site_index, other_species)) {
        const int target_particle_index = this->r.particle_index_at_position(target_site_index, other_species);
        BOOST_ASSERT(target_particle_index >= 0);
        move.push_back(SingleParticleMove(Particle(target_particle_index, other_species), this->r[particle]));
    }

    return move;
}

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) template class BCSWavefunction<type>
#include "vmc-supported-types.hpp"
