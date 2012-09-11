#include <boost/assert.hpp>
#include <boost/make_shared.hpp>
#include <boost/cast.hpp>

#include "vmc-math-utils.hpp"
#include "DBLWavefunction.hpp"

DBLWavefunction::Amplitude::Amplitude (const boost::shared_ptr<DBLWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction::Amplitude(wf_, r_),
      m_partial_update_step(0)
{
    BOOST_ASSERT(r.get_N_species() == 1);
    BOOST_ASSERT(r.get_N_sites() == wf_->orbital_def1->get_N_sites());
    BOOST_ASSERT(r.get_N_sites() == wf_->orbital_def2->get_N_sites());
    BOOST_ASSERT(r.get_N_filled(0) == wf_->orbital_def1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) == wf_->orbital_def2->get_N_filled());
    BOOST_ASSERT(wf_->orbital_def1->get_lattice_ptr() == wf_->orbital_def2->get_lattice_ptr());

    reinitialize();
}

void DBLWavefunction::Amplitude::perform_move_ (const Move &move)
{
    // we require that m_partial_update_step == 0 between moves; otherwise, psi_() will
    // return zero when it shouldn't.
    BOOST_ASSERT(m_partial_update_step == 0);

    m_current_move = move;

    do_perform_move<true>(move);
}

// During perform_move(), we want to short-circuit out as soon as either of the
// determinants is zero, as we already know what psi_() will be.  However,
// doing so means we must later make a second pass to finish the (possibly)
// remaining determinant update if finish_update() is called.  This templated
// function allows us to use the same code in both passes, with only minor
// differences.
template <bool first_pass>
void DBLWavefunction::Amplitude::do_perform_move (const Move &move)
{
    if (!first_pass && m_partial_update_step == 0)
        return;

    const DBLWavefunction *wf_ = boost::polymorphic_downcast<const DBLWavefunction *>(wf.get());

    switch (first_pass ? 2 : m_partial_update_step)
    {
    case 2:
        {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> d2_cols;
            for (unsigned int i = 0; i < move.size(); ++i)
                d2_cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
            cmat2.update_columns(d2_cols, wf_->orbital_def2->get_orbitals());
        }
        if (first_pass && cmat2.get_determinant() == amplitude_t(0)) {
            m_partial_update_step = 1;
            return;
        }

    case 1:
        {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> d1_cols;
            for (unsigned int i = 0; i < move.size(); ++i)
                d1_cols.push_back(std::make_pair(move[i].particle.index, move[i].destination));
            cmat1.update_columns(d1_cols, wf_->orbital_def1->get_orbitals());
        }
    }

    if (!first_pass)
        m_partial_update_step = 0;
}

amplitude_t DBLWavefunction::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return amplitude_t(0);

    // fixme: we could cache or precalculate this ... but i doubt it would make
    // much difference really
    const DBLWavefunction *wf_ = boost::polymorphic_downcast<const DBLWavefunction *>(wf.get());
    return (complex_pow(cmat1.get_determinant(), wf_->d1_exponent)
            * complex_pow(cmat2.get_determinant(), wf_->d2_exponent));
}

void DBLWavefunction::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    cmat1.finish_columns_update();
    cmat2.finish_columns_update();
}

void DBLWavefunction::Amplitude::cancel_move_ (void)
{
    switch (m_partial_update_step) {
    case 0:
        cmat1.cancel_columns_update();
    case 1:
        cmat2.cancel_columns_update();
    }

    m_partial_update_step = 0;
}

void DBLWavefunction::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    (void) species; // will always be 0
    cmat1.swap_columns(particle1_index, particle2_index);
    cmat2.swap_columns(particle1_index, particle2_index);
}

void DBLWavefunction::Amplitude::reset_ (const PositionArguments &r_)
{
    const DBLWavefunction *wf_ = boost::polymorphic_downcast<const DBLWavefunction *>(wf.get());

    BOOST_ASSERT(r_.get_N_species() == 1);
    BOOST_ASSERT(r_.get_N_sites() == wf_->orbital_def1->get_N_sites());
    BOOST_ASSERT(r_.get_N_filled(0) == wf_->orbital_def1->get_N_filled());

    r = r_;
    reinitialize();
}

void DBLWavefunction::Amplitude::reinitialize (void)
{
    const DBLWavefunction *wf_ = boost::polymorphic_downcast<const DBLWavefunction *>(wf.get());

    const unsigned int N = r.get_N_filled(0);
    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat1(N, N), mat2(N, N);
    for (unsigned int i = 0; i < N; ++i) {
        const Particle particle(i, 0);
        mat1.col(i) = wf_->orbital_def1->at_position(r[particle]);
        mat2.col(i) = wf_->orbital_def2->at_position(r[particle]);
    }
    cmat1 = mat1;
    cmat2 = mat2;
}

boost::shared_ptr<Wavefunction::Amplitude> DBLWavefunction::Amplitude::clone_ (void) const
{
    return boost::make_shared<DBLWavefunction::Amplitude>(*this);
}
