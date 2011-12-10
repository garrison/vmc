#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "DBMWavefunctionAmplitude.hpp"

DBMWavefunctionAmplitude::DBMWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_d1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_d2, const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_up, const boost::shared_ptr<const OrbitalDefinitions> &orbital_f_down)
    : WavefunctionAmplitude(r_, orbital_d1->get_lattice_ptr()),
      m_orbital_d1(orbital_d1),
      m_orbital_d2(orbital_d2),
      m_orbital_f_up(orbital_f_up),
      m_orbital_f_down(orbital_f_down)
{
    BOOST_ASSERT(orbital_d1->get_lattice_ptr() == orbital_d2->get_lattice_ptr());
    BOOST_ASSERT(orbital_d1->get_lattice_ptr() == orbital_f_up->get_lattice_ptr());
    BOOST_ASSERT(orbital_d1->get_lattice_ptr() == orbital_f_down->get_lattice_ptr());

    BOOST_ASSERT(r.get_N_sites() == orbital_d1->get_N_sites());

    BOOST_ASSERT(r.get_N_filled() == orbital_d1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled() == orbital_d2->get_N_filled());
    BOOST_ASSERT(r.get_N_filled() == orbital_f_up->get_N_filled() + orbital_f_down->get_N_filled());

    reinitialize();
}

void DBMWavefunctionAmplitude::move_particle_ (unsigned int particle, unsigned int new_site_index)
{
    unsigned int N = r.get_N_filled();
    unsigned int M = m_orbital_f_up->get_N_filled();

    BOOST_ASSERT(particle < N);

    r.update_position(particle, new_site_index);

    // update the Ceperley matrices
    m_particle_moved_is_up = bool(particle < M);
    m_cmat_d1.update_column(particle, m_orbital_d1->at_position(new_site_index));
    m_cmat_d2.update_column(particle, m_orbital_d2->at_position(new_site_index));
    if (m_particle_moved_is_up)
        m_cmat_f_up.update_column(particle, m_orbital_f_up->at_position(new_site_index));
    else
        m_cmat_f_down.update_column(particle - M, m_orbital_f_down->at_position(new_site_index));
}

amplitude_t DBMWavefunctionAmplitude::psi_ (void) const
{
    return (m_cmat_d1.get_determinant()
            * m_cmat_d2.get_determinant()
            * m_cmat_f_up.get_determinant()
            * m_cmat_f_down.get_determinant());
}

void DBMWavefunctionAmplitude::finish_particle_moved_update_ (void)
{
    m_cmat_d1.finish_column_update();
    m_cmat_d2.finish_column_update();
    (m_particle_moved_is_up ? m_cmat_f_up : m_cmat_f_down).finish_column_update();
}

void DBMWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    BOOST_ASSERT(r_.get_N_sites() == m_orbital_d1->get_N_sites());
    BOOST_ASSERT(r_.get_N_filled() == m_orbital_d1->get_N_filled());

    r = r_;
    reinitialize();
}

void DBMWavefunctionAmplitude::reinitialize (void)
{
    const unsigned int N = r.get_N_filled();
    const unsigned int M = m_orbital_f_up->get_N_filled();

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_d1(N, N);
    for (unsigned int i = 0; i < N; ++i)
        mat_d1.col(i) = m_orbital_d1->at_position(r[i]);
    m_cmat_d1 = mat_d1;

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_d2(N, N);
    for (unsigned int i = 0; i < N; ++i)
        mat_d2.col(i) = m_orbital_d2->at_position(r[i]);
    m_cmat_d2 = mat_d2;

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_f_up(M, M);
    for (unsigned int i = 0; i < M; ++i)
        mat_f_up.col(i) = m_orbital_f_up->at_position(r[i]);
    m_cmat_f_up = mat_f_up;

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat_f_down(N - M, N - M);
    for (unsigned int i = 0; i < N - M; ++i)
        mat_f_down.col(i) = m_orbital_f_down->at_position(r[i + M]);
    m_cmat_f_down = mat_f_down;
}

boost::shared_ptr<WavefunctionAmplitude> DBMWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<DBMWavefunctionAmplitude>(*this);
}
