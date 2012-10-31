#include <vector>
#include <algorithm>

#include <boost/make_shared.hpp>
#include <boost/cast.hpp>

#include "BCSWavefunction.hpp"
#include "RandomNumberGenerator.hpp"
#include "random-combination.hpp"
#include "random-configuration.hpp"
#include "random-move.hpp"

BCSWavefunction::Amplitude::Amplitude (const boost::shared_ptr<const BCSWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction::Amplitude(wf_, r_),
      m_current_jastrow(1),
      m_partial_update_step(0)
{    
    reinitialize();
}

void BCSWavefunction::Amplitude::perform_move_ (const Move &move)
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
template <bool first_pass>
void BCSWavefunction::Amplitude::do_perform_move (const Move &move)
{
    const BCSWavefunction *wf_ = boost::polymorphic_downcast<const BCSWavefunction *>(wf.get());
    const boost::shared_ptr<const Lattice> &lattice = wf->lattice;
    const unsigned int M = wf_->M;

    if (first_pass) {
        // explicitly enforce Gutzwiller projection
        for (unsigned int i = 0; i < move.size(); ++i) {
            if (r.is_occupied(move[i].destination, move[i].particle.species ^ 1)) {
                m_partial_update_step = 2;
                return;
            }
        }
    }

    switch (first_pass ? 2 : m_partial_update_step) {
    case 2:
        if (wf_->jastrow) {
            m_current_jastrow = wf_->jastrow->compute_jastrow(r);
            if (m_current_jastrow == amplitude_t(0)) {
                m_partial_update_step = 1;
                return;
            }
        } else {
            BOOST_ASSERT(m_current_jastrow == 1);
        }

    case 1:
        {
            Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> srcmat(M, M);
            lw_vector<unsigned int, MAX_MOVE_SIZE> rows, cols;

            for (unsigned int j = 0; j < move.size(); ++j) {
                const unsigned int particle_index = move[j].particle.index;
                const LatticeSite new_site(lattice->site_from_index(move[j].destination));

                if (move[j].particle.species == 0) {
                    const std::vector<unsigned int> & dn_pos = r.r_vector(1);
                    for (unsigned int i = 0; i < M; ++i) {
                        LatticeSite rup_minus_rdn(new_site);
                        lattice->asm_subtract_site_vector(rup_minus_rdn, lattice->site_from_index(dn_pos[i]).bravais_site());
                        lattice->enforce_boundary(rup_minus_rdn);
                        srcmat(particle_index, i) = wf_->phi[lattice->site_to_index(rup_minus_rdn)];
                    }
                    rows.push_back(particle_index);
                } else {
                    BOOST_ASSERT(move[j].particle.species == 1);
                    const std::vector<unsigned int> & up_pos = r.r_vector(0);
                    for (unsigned int i = 0; i < M; ++i) {
                        LatticeSite rup_minus_rdn(lattice->site_from_index(up_pos[i]));
                        lattice->asm_subtract_site_vector(rup_minus_rdn, new_site.bravais_site());
                        lattice->enforce_boundary(rup_minus_rdn);
                        srcmat(i, particle_index) = wf_->phi[lattice->site_to_index(rup_minus_rdn)];
                    }
                    cols.push_back(particle_index);
                }
            }

            m_cmat.update_rows_and_columns(rows, cols, srcmat);
        }

        m_partial_update_step = 0;

    case 0: ;
    }
}

amplitude_t BCSWavefunction::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return amplitude_t(0);

    return m_current_jastrow * m_cmat.get_determinant();
}

void BCSWavefunction::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    m_cmat.finish_rows_and_columns_update();
}

void BCSWavefunction::Amplitude::cancel_move_ (void)
{
    if (m_partial_update_step == 0)
        m_cmat.cancel_rows_and_columns_update();
    m_current_jastrow = m_old_jastrow; // this is always safe as a way to reset
                                       // m_current_jastrow, even during a
                                       // partial update
    m_partial_update_step = 0;
}

void BCSWavefunction::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    if (species == 0) {
        m_cmat.swap_rows(particle1_index, particle2_index);
    } else {
        BOOST_ASSERT(species == 1);
        m_cmat.swap_columns(particle1_index, particle2_index);
    }
}

void BCSWavefunction::Amplitude::reset_ (const PositionArguments &r_)
{
    r = r_;
    reinitialize();
}

void BCSWavefunction::Amplitude::reinitialize (void)
{
    const BCSWavefunction *wf_ = boost::polymorphic_downcast<const BCSWavefunction *>(wf.get());
    const boost::shared_ptr<const Lattice> &lattice = wf->lattice;
    const unsigned int M = wf_->M;

    BOOST_ASSERT(r.get_N_sites() == lattice->total_sites());
    BOOST_ASSERT(r.get_N_species() == 2);
    BOOST_ASSERT(r.get_N_filled(0) == M);
    BOOST_ASSERT(r.get_N_filled(1) == M);

    if (wf_->jastrow) {
        m_current_jastrow = wf_->jastrow->compute_jastrow(r);
    }

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_phi(M, M);
    const std::vector<unsigned int> & up_pos = r.r_vector(0);
    const std::vector<unsigned int> & dn_pos = r.r_vector(1);
    for (unsigned int i = 0; i < M; ++i) {
        const LatticeSite rup(lattice->site_from_index(up_pos[i]));
        for (unsigned int j = 0; j < M; ++j) {
            LatticeSite rup_minus_rdn(rup);
            lattice->asm_subtract_site_vector(rup_minus_rdn, lattice->site_from_index(dn_pos[j]).bravais_site());
            lattice->enforce_boundary(rup_minus_rdn);
            mat_phi(i, j) = wf_->phi[lattice->site_to_index(rup_minus_rdn)];
        }
    }

    m_cmat = mat_phi;
}

boost::shared_ptr<Wavefunction::Amplitude> BCSWavefunction::Amplitude::clone_ (void) const
{
    return boost::make_shared<BCSWavefunction::Amplitude>(*this);
}

boost::shared_ptr<Wavefunction::Amplitude> BCSWavefunction::create_nonzero_wavefunctionamplitude (const boost::shared_ptr<const Wavefunction> &this_ptr, RandomNumberGenerator &rng, unsigned int n_attempts) const
{
    const unsigned int M = get_N_filled(0);
    const unsigned int N = lattice->total_sites();

    BOOST_ASSERT(2 * M <= N);

    while (n_attempts--) {
        // Take into account Gutzwiller projection [at most one spinon ("particle") per site]
        std::vector<std::vector<unsigned int> > vv(2);
        vv[0] = some_random_configuration(M, *lattice, rng);
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

        boost::shared_ptr<Wavefunction::Amplitude> wfa(create_wavefunctionamplitude(this_ptr, PositionArguments(vv, lattice->total_sites())));
        if (wfa->psi() != amplitude_t(0))
            return wfa;
    }
    return boost::shared_ptr<Wavefunction::Amplitude>();
}

Move BCSWavefunction::Amplitude::propose_move (RandomNumberGenerator &rng) const
{
    const BCSWavefunction *wf_ = boost::polymorphic_downcast<const BCSWavefunction *>(wf.get());
    const unsigned int N = wf->lattice->total_sites();
    const unsigned int M = wf_->M;

    BOOST_ASSERT(2 * M <= N);
    const unsigned int empty_sites = N - 2 * M;

    if (empty_sites != 0) {
        if (rng.random_small_uint(2) == 0) { // FIXME: consider a more efficient heuristic
            // move a random particle to an empty site (by calling the
            // superclass's implementation)
            return Wavefunction::Amplitude::propose_move(rng);
        }
    }

    // choose a particle of each species and swap their positions
    Move move;
    const Particle up_particle(rng.random_small_uint(M), 0);
    const Particle dn_particle(rng.random_small_uint(M), 1);
    move.push_back(SingleParticleMove(up_particle, r[dn_particle]));
    move.push_back(SingleParticleMove(dn_particle, r[up_particle]));
    return move;
}
