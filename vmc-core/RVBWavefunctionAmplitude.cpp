#include <vector>
#include <algorithm>

#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "RVBWavefunctionAmplitude.hpp"

RVBWavefunctionAmplitude::RVBWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const Lattice> &lattice_, const std::vector<complex_t> &phi)
    : WavefunctionAmplitude(r_, lattice_),
      m_update_in_progress(false),
      m_phi(phi)
{    
    reinitialize();
}

void RVBWavefunctionAmplitude::perform_move_ (Particle particle, unsigned int new_site_index)
{
    BOOST_ASSERT(r.particle_is_valid(particle));
    BOOST_ASSERT(new_site_index < r.get_N_sites());

    const unsigned int species = particle.species;
    const unsigned int old_site_index = r[particle];

    r.update_position(particle, new_site_index);

    const unsigned int other_species = (species == 0) ? 1 : 0;
    const int other_particle_index = r.particle_index_at_pos(new_site_index, other_species);
    BOOST_ASSERT(other_particle_index != -1);
    const Particle other_particle(other_particle_index, other_species);

    r.update_position(other_particle, old_site_index);

    const unsigned int M = r.get_N_filled(0);
    BOOST_ASSERT(M == r.get_N_filled(1));

    const unsigned int moved_up_particle_index = (particle.species == 0) ? particle.index : other_particle.index;
    const unsigned int moved_down_particle_index = (particle.species == 1) ? particle.index : other_particle.index;

    const LatticeSite new_site_for_up(lattice->site_from_index(particle.species == 0 ? new_site_index : old_site_index));
    const LatticeSite new_site_for_down(lattice->site_from_index(particle.species == 1 ? new_site_index : old_site_index));

    Eigen::Matrix<complex_t, Eigen::Dynamic, 1> new_row(M);
    const std::vector<unsigned int> & down_pos = r.r_vector(1);
    for (unsigned int i = 0; i < M; ++i) {
        LatticeSite rup_minus_rdown(new_site_for_up);
        lattice->asm_subtract_site_vector(rup_minus_rdown, lattice->site_from_index(down_pos[i]).bravais_site());
        new_row[i] = m_phi[lattice->site_to_index(rup_minus_rdown)];
    }

    m_update_in_progress = true;
    m_new_cmat = m_cmat;

    m_new_cmat.update_row(moved_up_particle_index, new_row);
    m_new_cmat.finish_row_update();
    
    Eigen::Matrix<complex_t, Eigen::Dynamic, 1> new_col(M);
    const std::vector<unsigned int> & up_pos = r.r_vector(0);
    for (unsigned int i = 0; i < M; ++i) {
        LatticeSite rup_minus_rdown(lattice->site_from_index(up_pos[i]));
        lattice->asm_subtract_site_vector(rup_minus_rdown, new_site_for_down.bravais_site());
        new_col[i] = m_phi[lattice->site_to_index(rup_minus_rdown)];
    }

    m_new_cmat.update_column(moved_down_particle_index, new_col);
}

amplitude_t RVBWavefunctionAmplitude::psi_ (void) const
{
    return (m_update_in_progress ? m_new_cmat : m_cmat).get_determinant();
}

void RVBWavefunctionAmplitude::finish_move_ (void)
{
    m_new_cmat.finish_column_update();
    m_cmat = m_new_cmat;
    m_update_in_progress = false;
}

void RVBWavefunctionAmplitude::cancel_move_ (Particle particle, unsigned int old_site_index)
{
    unsigned int other_species = particle.species ^ 1;
    const int other_particle_index = r.particle_index_at_pos(old_site_index, other_species);
    BOOST_ASSERT(other_particle_index >= 0);
    r.update_position(Particle(other_particle_index, other_species), r[particle]);
    r.update_position(particle, old_site_index);

    m_update_in_progress = false;
}

void RVBWavefunctionAmplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    // we don't support Renyi calculations in RVBWavefunctionAmplitude yet
    BOOST_ASSERT(false);

    (void) particle1_index;
    (void) particle2_index;
    (void) species;
}

void RVBWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    r = r_;
    reinitialize();
}

void RVBWavefunctionAmplitude::reinitialize (void)
{
    BOOST_ASSERT(r.get_N_species() == 2);

    BOOST_ASSERT(r.get_N_sites() == lattice->total_sites());

    // Each species is at half-filling and assume unpolarized wave function
    BOOST_ASSERT(2 * r.get_N_filled(0) == lattice->total_sites());
    BOOST_ASSERT(2 * r.get_N_filled(1) == lattice->total_sites());

    BOOST_ASSERT(r.get_N_sites() == m_phi.size());

    const unsigned int M = r.get_N_filled(0);

    Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic> mat_phi(M, M);
    const std::vector<unsigned int> & up_pos = r.r_vector(0);
    const std::vector<unsigned int> & down_pos = r.r_vector(1);
    for (unsigned int i = 0; i < M; ++i) {
        const LatticeSite rup(lattice->site_from_index(up_pos[i]));
        for (unsigned int j = 0; j < M; ++j) {
            LatticeSite rup_minus_rdown(rup);
            lattice->asm_subtract_site_vector(rup_minus_rdown, lattice->site_from_index(down_pos[j]).bravais_site());
            mat_phi(i, j) = m_phi[lattice->site_to_index(rup_minus_rdown)];
        }
    }

    m_cmat = mat_phi;
}

boost::shared_ptr<WavefunctionAmplitude> RVBWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<RVBWavefunctionAmplitude>(*this);
}

void RVBWavefunctionAmplitude::reset_with_filler (const RandomFiller &filler, rng_class &rng)
{
    BOOST_ASSERT(r.get_N_species() == 2);

    const unsigned int M = r.get_N_filled(0);
    BOOST_ASSERT(M == r.get_N_filled(1));

    const unsigned int N_sites = lattice->total_sites();
    BOOST_ASSERT(N_sites == r.get_N_sites());

    BOOST_ASSERT(r.get_N_filled_total() == N_sites);  // Assert a spin wave function!

    // Take into account Gutzwiller projection [exactly one spinon ("particle") per site]
    std::vector<std::vector<unsigned int> > vv(2);
    vv[0] = filler.some_random_filling(M, rng);
    // NOTE: this method requires O(N ^ 2) time, but this could technically be
    // done in O(N) time.  Not a big deal here.
    for (unsigned int i = 0; i < N_sites; ++i) {
        std::vector<unsigned int>::iterator it = std::find(vv[0].begin(), vv[0].end(), i);
        if (it == vv[0].end())
            vv[1].push_back(i);
    }

    BOOST_ASSERT(vv[0].size() == vv[1].size());
    reset(PositionArguments(vv, N_sites));
}
