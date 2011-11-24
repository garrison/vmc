#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "FreeFermionWavefunctionAmplitude.hpp"

FreeFermionWavefunctionAmplitude::FreeFermionWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_)
    : WavefunctionAmplitude(r_, orbital_def_->get_lattice_ptr()),
      orbital_def(orbital_def_)
{
    BOOST_ASSERT(r.get_N_sites() == orbital_def->get_N_sites());
    BOOST_ASSERT(r.get_N_filled() == orbital_def->get_N_filled());

    reinitialize();
}

void FreeFermionWavefunctionAmplitude::move_particle_ (unsigned int particle, unsigned int new_site_index)
{
    unsigned int N = r.get_N_filled();
    BOOST_ASSERT(particle < N);

    r.update_position(particle, new_site_index);

    // calculate each phi at new position and update the Ceperley matrix
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, 1> phivec(N);
    cmat.update_row(particle, orbital_def->get(new_site_index));
}

amplitude_t FreeFermionWavefunctionAmplitude::psi_ (void) const
{
    return cmat.get_determinant();
}

void FreeFermionWavefunctionAmplitude::finish_particle_moved_update_ (void)
{
    cmat.finish_row_update();
}

void FreeFermionWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    BOOST_ASSERT(r.get_N_sites() == orbital_def->get_N_sites());
    BOOST_ASSERT(r.get_N_filled() == orbital_def->get_N_filled());

    r = r_;
    reinitialize();
}

void FreeFermionWavefunctionAmplitude::reinitialize (void)
{
    unsigned int N = r.get_N_filled();
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    for (unsigned int i = 0; i < N; ++i)
        mat.row(i) = orbital_def->get(r[i]);
    cmat = mat;
}

boost::shared_ptr<WavefunctionAmplitude> FreeFermionWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<FreeFermionWavefunctionAmplitude>(*this);
}
