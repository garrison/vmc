#include <vector>

#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "vmc-math-utils.hpp"
#include "DMetalWavefunctionAmplitude.hpp"
#include "random-filling.hpp"

DMetalWavefunctionAmplitude::DMetalWavefunctionAmplitude (const PositionArguments &r_,
                                                          const boost::shared_ptr<const OrbitalDefinitions> &orbital_d1,
                                                          const boost::shared_ptr<const OrbitalDefinitions> &orbital_d2,
                                                          const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_up,
                                                          const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_down,
                                                          real_t d1_exponent,
                                                          real_t d2_exponent,
                                                          real_t f_up_exponent,
                                                          real_t f_down_exponent)
    : WavefunctionAmplitude(r_, orbital_d1->get_lattice_ptr()),
      m_orbital_d1(orbital_d1),
      m_orbital_d2(orbital_d2),
      m_orbital_f_up(orbital_f_up),
      m_orbital_f_down(orbital_f_down),
      m_d1_exponent(d1_exponent),
      m_d2_exponent(d2_exponent),
      m_f_up_exponent(f_up_exponent),
      m_f_down_exponent(f_down_exponent),
      m_gutzwiller_rejection_in_progress(false)
{
    reinitialize();
}

void DMetalWavefunctionAmplitude::perform_move_ (Particle particle, unsigned int new_site_index)
{
    BOOST_ASSERT(r.particle_is_valid(particle));
    BOOST_ASSERT(new_site_index < r.get_N_sites());

    r.update_position(particle, new_site_index);

    BOOST_ASSERT(!m_gutzwiller_rejection_in_progress);
    if (r.is_occupied(new_site_index, particle.species ^ 1)) {
        // explicitly enforce Gutzwiller projection
        m_gutzwiller_rejection_in_progress = true;
        return;
    }

    const unsigned int M = m_orbital_f_up->get_N_filled();

    // update the Ceperley matrices
    m_particle_moved_is_up = bool(particle.species == 0);
    const unsigned int particle_column_index = m_particle_moved_is_up ? particle.index : particle.index + M;
    m_cmat_d1.update_column(particle_column_index, m_orbital_d1->at_position(new_site_index));
    m_cmat_d2.update_column(particle_column_index, m_orbital_d2->at_position(new_site_index));
    if (m_particle_moved_is_up)
        m_cmat_f_up.update_column(particle.index, m_orbital_f_up->at_position(new_site_index));
    else
        m_cmat_f_down.update_column(particle.index, m_orbital_f_down->at_position(new_site_index));
}

amplitude_t DMetalWavefunctionAmplitude::psi_ (void) const
{
    if (m_gutzwiller_rejection_in_progress)
        return amplitude_t(0);

    return (complex_pow(m_cmat_d1.get_determinant(), m_d1_exponent)
            * complex_pow(m_cmat_d2.get_determinant(), m_d2_exponent)
            * complex_pow(m_cmat_f_up.get_determinant(), m_f_up_exponent)
            * complex_pow(m_cmat_f_down.get_determinant(), m_f_down_exponent));
}

void DMetalWavefunctionAmplitude::finish_move_ (void)
{
    if (m_gutzwiller_rejection_in_progress) {
        m_gutzwiller_rejection_in_progress = false;
        return;
    }

    m_cmat_d1.finish_column_update();
    m_cmat_d2.finish_column_update();
    (m_particle_moved_is_up ? m_cmat_f_up : m_cmat_f_down).finish_column_update();
}

void DMetalWavefunctionAmplitude::cancel_move_ (Particle particle, unsigned int old_site_index)
{
    r.update_position(particle, old_site_index);

    if (m_gutzwiller_rejection_in_progress) {
        m_gutzwiller_rejection_in_progress = false;
        return;
    }

    m_cmat_d1.cancel_column_update();
    m_cmat_d2.cancel_column_update();
    (m_particle_moved_is_up ? m_cmat_f_up : m_cmat_f_down).cancel_column_update();
}

void DMetalWavefunctionAmplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    const unsigned int M = m_orbital_f_up->get_N_filled();
    const unsigned int particle1_column_index = (species == 0) ? particle1_index : particle1_index + M;
    const unsigned int particle2_column_index = (species == 0) ? particle2_index : particle2_index + M;
    m_cmat_d1.swap_columns(particle1_column_index, particle2_column_index);
    m_cmat_d2.swap_columns(particle1_column_index, particle2_column_index);
    (species == 0 ? m_cmat_f_up : m_cmat_f_down).swap_columns(particle1_index, particle2_index);
}

void DMetalWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    r = r_;
    reinitialize();
}

void DMetalWavefunctionAmplitude::reinitialize (void)
{
    BOOST_ASSERT(r.get_N_species() == 2);

    BOOST_ASSERT(r.get_N_sites() == m_orbital_d1->get_N_sites());
    BOOST_ASSERT(m_orbital_d1->get_lattice_ptr() == m_orbital_d2->get_lattice_ptr());
    BOOST_ASSERT(m_orbital_d1->get_lattice_ptr() == m_orbital_f_up->get_lattice_ptr());
    BOOST_ASSERT(m_orbital_d1->get_lattice_ptr() == m_orbital_f_down->get_lattice_ptr());

    BOOST_ASSERT(r.get_N_filled(0) + r.get_N_filled(1) == m_orbital_d1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) + r.get_N_filled(1) == m_orbital_d2->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) == m_orbital_f_up->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(1) == m_orbital_f_down->get_N_filled());

    const unsigned int N = m_orbital_d1->get_N_filled();
    const unsigned int M = m_orbital_f_up->get_N_filled();

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_d1(N, N);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_d2(N, N);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_f_up(M, M);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_f_down(N - M, N - M);

    for (unsigned int i = 0; i < r.get_N_filled(0); ++i) {
        const Particle particle(i, 0);
        mat_d1.col(i) = m_orbital_d1->at_position(r[particle]);
        mat_d2.col(i) = m_orbital_d2->at_position(r[particle]);
        mat_f_up.col(i) = m_orbital_f_up->at_position(r[particle]);
    }

    for (unsigned int i = 0; i < r.get_N_filled(1); ++i) {
        const Particle particle(i, 1);
        mat_d1.col(i + M) = m_orbital_d1->at_position(r[particle]);
        mat_d2.col(i + M) = m_orbital_d2->at_position(r[particle]);
        mat_f_down.col(i) = m_orbital_f_down->at_position(r[particle]);
    }

    m_cmat_d1 = mat_d1;
    m_cmat_d2 = mat_d2;
    m_cmat_f_up = mat_f_up;
    m_cmat_f_down = mat_f_down;
}

boost::shared_ptr<WavefunctionAmplitude> DMetalWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<DMetalWavefunctionAmplitude>(*this);
}

void DMetalWavefunctionAmplitude::reset_with_random_positions (rng_class &rng)
{
    const unsigned int N = m_orbital_d1->get_N_filled();
    const unsigned int M = m_orbital_f_up->get_N_filled();

    // take into account the Gutzwiller projection
    std::vector<std::vector<unsigned int> > vv(2);
    vv[0] = some_random_filling(N, *lattice, rng);
    for (unsigned int i = M; i < N; ++i)
        vv[1].push_back(vv[0][i]);
    vv[0].resize(M);

    reset(PositionArguments(vv, lattice->total_sites()));
}
