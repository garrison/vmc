#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "FreeFermionWavefunctionAmplitude.hpp"
#include "HypercubicLattice.hpp"

FreeFermionWavefunctionAmplitude::FreeFermionWavefunctionAmplitude (const PositionArguments &r_)
    : WavefunctionAmplitude(r_),
      orbital_def(r_.get_N_filled())
{
    const boost::array<int, 1> a = { { r_.get_N_sites() } };
    lattice = boost::make_shared<const HypercubicLattice<1> >(a);
    reinitialize();
}

void FreeFermionWavefunctionAmplitude::move_particle_ (unsigned int particle, unsigned int new_site_index)
{
    unsigned int N = r.get_N_filled();
    BOOST_ASSERT(particle < N);

    r.update_position(particle, new_site_index);

    // calculate each phi at new position and update the Ceperley matrix
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, 1> phivec(N);
    for (unsigned int i = 0; i < N; ++i)
	phivec(i) = orbital_def.phi(i, new_site_index, *lattice);
    cmat.update_row(particle, phivec);
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
    r = r_;
    reinitialize();
}

void FreeFermionWavefunctionAmplitude::reinitialize (void)
{
    unsigned int N = r.get_N_filled();
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    for (unsigned int i = 0; i < N; ++i) {
	for (unsigned int j = 0; j < N; ++j)
	    mat(i, j) = orbital_def.phi(j, r[i], *lattice);
    }
    cmat = mat;
}

boost::shared_ptr<WavefunctionAmplitude> FreeFermionWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<FreeFermionWavefunctionAmplitude>(*this);
}
