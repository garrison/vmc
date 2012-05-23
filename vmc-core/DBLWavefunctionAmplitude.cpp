#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "vmc-math-utils.hpp"
#include "DBLWavefunctionAmplitude.hpp"

DBLWavefunctionAmplitude::DBLWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_2, real_t d1_exponent_, real_t d2_exponent_)
    : WavefunctionAmplitude(r_, orbital_def_1->get_lattice_ptr()),
      orbital_def1(orbital_def_1),
      orbital_def2(orbital_def_2),
      d1_exponent(d1_exponent_),
      d2_exponent(d2_exponent_)
{
    BOOST_ASSERT(r.get_N_species() == 1);
    BOOST_ASSERT(r.get_N_sites() == orbital_def1->get_N_sites());
    BOOST_ASSERT(r.get_N_sites() == orbital_def2->get_N_sites());
    BOOST_ASSERT(r.get_N_filled(0) == orbital_def1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) == orbital_def2->get_N_filled());
    BOOST_ASSERT(orbital_def1->get_lattice_ptr() == orbital_def2->get_lattice_ptr());

    reinitialize();
}

void DBLWavefunctionAmplitude::perform_move_ (Particle particle, unsigned int new_site_index)
{
    BOOST_ASSERT(r.particle_is_valid(particle));
    BOOST_ASSERT(new_site_index < r.get_N_sites());

    r.update_position(particle, new_site_index);

    // update the Ceperley matrices
    cmat1.update_column(particle.index, orbital_def1->at_position(new_site_index));
    cmat2.update_column(particle.index, orbital_def2->at_position(new_site_index));
}

amplitude_t DBLWavefunctionAmplitude::psi_ (void) const
{
    // fixme: we could cache or precalculate this ... but i doubt it would make
    // much difference really
    return (complex_pow(cmat1.get_determinant(), d1_exponent)
            * complex_pow(cmat2.get_determinant(), d2_exponent));
}

void DBLWavefunctionAmplitude::finish_move_ (void)
{
    cmat1.finish_column_update();
    cmat2.finish_column_update();
}

void DBLWavefunctionAmplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    (void) species; // will always be 0
    cmat1.swap_columns(particle1_index, particle2_index);
    cmat2.swap_columns(particle1_index, particle2_index);
}

void DBLWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    BOOST_ASSERT(r_.get_N_species() == 1);
    BOOST_ASSERT(r_.get_N_sites() == orbital_def1->get_N_sites());
    BOOST_ASSERT(r_.get_N_filled(0) == orbital_def1->get_N_filled());

    r = r_;
    reinitialize();
}

void DBLWavefunctionAmplitude::reinitialize (void)
{
    const unsigned int N = r.get_N_filled(0);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat1(N, N), mat2(N, N);
    for (unsigned int i = 0; i < N; ++i) {
        const Particle particle(i, 0);
        mat1.col(i) = orbital_def1->at_position(r[particle]);
        mat2.col(i) = orbital_def2->at_position(r[particle]);
    }
    cmat1 = mat1;
    cmat2 = mat2;
}

boost::shared_ptr<WavefunctionAmplitude> DBLWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<DBLWavefunctionAmplitude>(*this);
}
