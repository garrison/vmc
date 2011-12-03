#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "DBLWavefunctionAmplitude.hpp"

DBLWavefunctionAmplitude::DBLWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_2)
    : WavefunctionAmplitude(r_, orbital_def_1->get_lattice_ptr()),
      orbital_def1(orbital_def_1),
      orbital_def2(orbital_def_2)
{
    BOOST_ASSERT(r.get_N_sites() == orbital_def1->get_N_sites());
    BOOST_ASSERT(r.get_N_sites() == orbital_def2->get_N_sites());
    BOOST_ASSERT(r.get_N_filled() == orbital_def1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled() == orbital_def2->get_N_filled());
    BOOST_ASSERT(orbital_def1->get_lattice_ptr() == orbital_def2->get_lattice_ptr());

    reinitialize();
}

void DBLWavefunctionAmplitude::move_particle_ (unsigned int particle, unsigned int new_site_index)
{
    unsigned int N = r.get_N_filled();
    BOOST_ASSERT(particle < N);

    r.update_position(particle, new_site_index);

    // update the Ceperley matrices
    cmat1.update_row(particle, orbital_def1->at_position(new_site_index));
    cmat2.update_row(particle, orbital_def2->at_position(new_site_index));
}

amplitude_t DBLWavefunctionAmplitude::psi_ (void) const
{
    return cmat1.get_determinant() * cmat2.get_determinant();
}

void DBLWavefunctionAmplitude::finish_particle_moved_update_ (void)
{
    cmat1.finish_row_update();
    cmat2.finish_row_update();
}

void DBLWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    BOOST_ASSERT(r_.get_N_sites() == orbital_def1->get_N_sites());
    BOOST_ASSERT(r_.get_N_filled() == orbital_def1->get_N_filled());

    r = r_;
    reinitialize();
}

void DBLWavefunctionAmplitude::reinitialize (void)
{
    unsigned int N = r.get_N_filled();
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat1(N, N), mat2(N, N);
    for (unsigned int i = 0; i < N; ++i) {
        mat1.row(i) = orbital_def1->at_position(r[i]);
        mat2.row(i) = orbital_def2->at_position(r[i]);
    }
    cmat1 = mat1;
    cmat2 = mat2;
}

boost::shared_ptr<WavefunctionAmplitude> DBLWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<DBLWavefunctionAmplitude>(*this);
}
