#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "FreeFermionWavefunctionAmplitude.hpp"

FreeFermionWavefunctionAmplitude::FreeFermionWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_)
    : WavefunctionAmplitude(r_, orbital_def_->get_lattice_ptr()),
      orbital_def(orbital_def_)
{
    BOOST_ASSERT(r.get_N_species() == 1);
    BOOST_ASSERT(r.get_N_sites() == orbital_def->get_N_sites());
    BOOST_ASSERT(r.get_N_filled(0) == orbital_def->get_N_filled());

    reinitialize();
}

void FreeFermionWavefunctionAmplitude::perform_move_ (const Move &move)
{
    lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> cols;
    for (unsigned int i = 0; i < move.size(); ++i)
        cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
    cmat.update_columns(cols, orbital_def->get_orbitals());
}

amplitude_t FreeFermionWavefunctionAmplitude::psi_ (void) const
{
    return cmat.get_determinant();
}

void FreeFermionWavefunctionAmplitude::finish_move_ (void)
{
    cmat.finish_columns_update();
}

void FreeFermionWavefunctionAmplitude::cancel_move_ (void)
{
    cmat.cancel_columns_update();
}

void FreeFermionWavefunctionAmplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    (void) species; // will always be zero
    cmat.swap_columns(particle1_index, particle2_index);
}

void FreeFermionWavefunctionAmplitude::reset_ (const PositionArguments &r_)
{
    BOOST_ASSERT(r_.get_N_species() == 1);
    BOOST_ASSERT(r_.get_N_sites() == orbital_def->get_N_sites());
    BOOST_ASSERT(r_.get_N_filled(0) == orbital_def->get_N_filled());

    r = r_;
    reinitialize();
}

void FreeFermionWavefunctionAmplitude::reinitialize (void)
{
    const unsigned int N = r.get_N_filled(0);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    for (unsigned int i = 0; i < N; ++i)
        mat.col(i) = orbital_def->at_position(r[Particle(i, 0)]);
    cmat = mat;
}

boost::shared_ptr<WavefunctionAmplitude> FreeFermionWavefunctionAmplitude::clone_ (void) const
{
    return boost::make_shared<FreeFermionWavefunctionAmplitude>(*this);
}
