#include <vector>
#include <algorithm>

#include <boost/assert.hpp>
#include <boost/make_shared.hpp>
#include <boost/cast.hpp>

#include "BCSWavefunction.hpp"
#include "random-configuration.hpp"
#include "random-move.hpp"

BCSWavefunction::Amplitude::Amplitude (const boost::shared_ptr<const BCSWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction::Amplitude(wf_, r_)
{    
    reinitialize();
}

void BCSWavefunction::Amplitude::perform_move_ (const Move &move)
{
    const BCSWavefunction *wf_ = boost::polymorphic_downcast<const BCSWavefunction *>(wf.get());
    const boost::shared_ptr<const Lattice> &lattice = wf->lattice;

    // first assert that it's a swap
    BOOST_ASSERT(move.size() == 2);
    BOOST_ASSERT(move[0].particle.species != move[1].particle.species);

    const unsigned int M = r.get_N_filled(0);
    BOOST_ASSERT(M == r.get_N_filled(1));

    // the next four lines are a bit opaque, but correct
    const unsigned int moved_up_particle_index = move[move[0].particle.species].particle.index;
    const unsigned int moved_dn_particle_index = move[move[1].particle.species].particle.index;
    const LatticeSite new_site_for_up(lattice->site_from_index(move[move[0].particle.species].destination));
    const LatticeSite new_site_for_dn(lattice->site_from_index(move[move[1].particle.species].destination));

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> srcmat(M, M);

    const std::vector<unsigned int> & dn_pos = r.r_vector(1);
    for (unsigned int i = 0; i < M; ++i) {
        LatticeSite rup_minus_rdn(new_site_for_up);
        lattice->asm_subtract_site_vector(rup_minus_rdn, lattice->site_from_index(dn_pos[i]).bravais_site());
        lattice->enforce_boundary(rup_minus_rdn);
        srcmat(moved_up_particle_index, i) = wf_->phi[lattice->site_to_index(rup_minus_rdn)];
    }

    const std::vector<unsigned int> & up_pos = r.r_vector(0);
    for (unsigned int i = 0; i < M; ++i) {
        LatticeSite rup_minus_rdn(lattice->site_from_index(up_pos[i]));
        lattice->asm_subtract_site_vector(rup_minus_rdn, new_site_for_dn.bravais_site());
        lattice->enforce_boundary(rup_minus_rdn);
        srcmat(i, moved_dn_particle_index) = wf_->phi[lattice->site_to_index(rup_minus_rdn)];
    }

    lw_vector<unsigned int, MAX_MOVE_SIZE> rows, cols;
    rows.push_back(moved_up_particle_index);
    cols.push_back(moved_dn_particle_index);
    m_cmat.update_rows_and_columns(rows, cols, srcmat);
}

amplitude_t BCSWavefunction::Amplitude::psi_ (void) const
{
    return m_cmat.get_determinant();
}

void BCSWavefunction::Amplitude::finish_move_ (void)
{
    m_cmat.finish_rows_and_columns_update();
}

void BCSWavefunction::Amplitude::cancel_move_ (void)
{
    m_cmat.cancel_rows_and_columns_update();
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

    BOOST_ASSERT(r.get_N_species() == 2);

    BOOST_ASSERT(r.get_N_sites() == lattice->total_sites());

    // Each species is at half-filling and assume unpolarized wave function
    BOOST_ASSERT(2 * r.get_N_filled(0) == lattice->total_sites());
    BOOST_ASSERT(2 * r.get_N_filled(1) == lattice->total_sites());

    BOOST_ASSERT(r.get_N_sites() == wf_->phi.size());

    const unsigned int M = r.get_N_filled(0);

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
    BOOST_ASSERT(2 * M == N);  // Assert a spin wave function!

    while (n_attempts--) {
        // Take into account Gutzwiller projection [exactly one spinon ("particle") per site]
        std::vector<std::vector<unsigned int> > vv(2);
        vv[0] = some_random_configuration(M, *lattice, rng);
        std::vector<unsigned int> occupied_sites(N);
        for (unsigned int i = 0; i < vv[0].size(); ++i) {
            BOOST_ASSERT(vv[0][i] < N);
            BOOST_ASSERT(occupied_sites[vv[0][i]] == 0);
            ++occupied_sites[vv[0][i]];
        }
        for (unsigned int i = 0; i < N; ++i) {
            if (occupied_sites[i] == 0)
                vv[1].push_back(i);
        }
        BOOST_ASSERT(vv[0].size() == vv[1].size());

        boost::shared_ptr<Wavefunction::Amplitude> wfa(create_wavefunctionamplitude(this_ptr, PositionArguments(vv, lattice->total_sites())));
        if (wfa->psi() != amplitude_t(0))
            return wfa;
    }
    return boost::shared_ptr<Wavefunction::Amplitude>();
}

Move BCSWavefunction::Amplitude::propose_move (RandomNumberGenerator &rng) const
{
    // choose a particle of each species and swap their positions
    Move move;
    const Particle particle(choose_random_particle(r, rng));
    const unsigned int other_species = particle.species ^ 1;
    const unsigned int proposed_site_index = plan_particle_move_to_nearby_empty_site(particle, r, *wf->lattice, rng);
    if (r.is_occupied(proposed_site_index, other_species)) {
        const int other_particle_index = r.particle_index_at_position(proposed_site_index, other_species);
        BOOST_ASSERT(other_particle_index >= 0);
        const Particle other_particle(other_particle_index, other_species);
        BOOST_ASSERT(r[other_particle] == proposed_site_index);

        move.push_back(SingleParticleMove(particle, r[other_particle]));
        move.push_back(SingleParticleMove(other_particle, r[particle]));
    }
    return move;
}
