#include <boost/assert.hpp>
#include <boost/make_shared.hpp>
#include <boost/cast.hpp>

#include "FreeFermionWavefunction.hpp"

FreeFermionWavefunction::Amplitude::Amplitude (const boost::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction::Amplitude(wf_, r_)
{
    BOOST_ASSERT(r.get_N_species() == 1);
    BOOST_ASSERT(r.get_N_sites() == wf_->orbital_def->get_N_sites());
    BOOST_ASSERT(r.get_N_filled(0) == wf_->orbital_def->get_N_filled());

    reinitialize();
}

void FreeFermionWavefunction::Amplitude::perform_move_ (const Move &move)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction *>(wf.get());

    lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> cols;
    for (unsigned int i = 0; i < move.size(); ++i)
        cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
    cmat.update_columns(cols, wf_->orbital_def->get_orbitals());
}

amplitude_t FreeFermionWavefunction::Amplitude::psi_ (void) const
{
    return cmat.get_determinant();
}

void FreeFermionWavefunction::Amplitude::finish_move_ (void)
{
    cmat.finish_columns_update();
}

void FreeFermionWavefunction::Amplitude::cancel_move_ (void)
{
    cmat.cancel_columns_update();
}

void FreeFermionWavefunction::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    (void) species; // will always be zero
    cmat.swap_columns(particle1_index, particle2_index);
}

void FreeFermionWavefunction::Amplitude::reset_ (const PositionArguments &r_)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction *>(wf.get());

    BOOST_ASSERT(r_.get_N_species() == 1);
    BOOST_ASSERT(r_.get_N_sites() == wf_->orbital_def->get_N_sites());
    BOOST_ASSERT(r_.get_N_filled(0) == wf_->orbital_def->get_N_filled());

    r = r_;
    reinitialize();
}

void FreeFermionWavefunction::Amplitude::reinitialize (void)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction *>(wf.get());

    const unsigned int N = r.get_N_filled(0);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
    for (unsigned int i = 0; i < N; ++i)
        mat.col(i) = wf_->orbital_def->at_position(r[Particle(i, 0)]);
    cmat = mat;
}

boost::shared_ptr<Wavefunction::Amplitude> FreeFermionWavefunction::Amplitude::clone_ (void) const
{
    return boost::make_shared<FreeFermionWavefunction::Amplitude>(*this);
}
