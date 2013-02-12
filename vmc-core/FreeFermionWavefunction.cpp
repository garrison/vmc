#include <boost/assert.hpp>
#include <boost/make_shared.hpp>
#include <boost/cast.hpp>

#include "FreeFermionWavefunction.hpp"

FreeFermionWavefunction::FreeFermionWavefunction (const std::vector<boost::shared_ptr<const OrbitalDefinitions> > &orbital_def_, const boost::shared_ptr<const JastrowFactor> &jastrow_)
    : Wavefunction(orbital_def_[0]->get_lattice_ptr()),
      orbital_def(orbital_def_),
      jastrow(jastrow_)
{
    BOOST_ASSERT(orbital_def.size() > 0);
#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    for (unsigned int i = 1; i < get_N_species(); ++i) {
        BOOST_ASSERT(orbital_def[0]->get_N_sites() == orbital_def[i]->get_N_sites());
    }
#endif
}

FreeFermionWavefunction::Amplitude::Amplitude (const boost::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction::Amplitude(wf_, r_),
      m_current_jastrow(1),
      m_partial_update_step(0),
      m_species_move_in_progress(wf_->get_N_species())
{
    reinitialize();
}

void FreeFermionWavefunction::Amplitude::perform_move_ (const Move &move)
{
    // we require that m_partial_update_step == 0 between moves; otherwise,
    // psi_() will return zero when it shouldn't.
    BOOST_ASSERT(m_partial_update_step == 0);

    // determine which species are being moved
    for (unsigned int i = 0; i < get_N_species(); ++i)
        m_species_move_in_progress[i] = false;
    for (unsigned int i = 0; i < move.size(); ++i)
        m_species_move_in_progress[move[i].particle.species] = true;

    m_current_move = move;

    m_old_jastrow = m_current_jastrow;

    do_perform_move<true>(move);
}

// During perform_move(), we want to short-circuit out as soon as any of the
// determinants is zero, as we already know what psi_() will be.  However,
// doing so means we must later make a second pass to finish all the
// determinant updates if finish_update() is called.  This templated function
// allows us to use the same code in both passes, with only minor differences.
template <bool first_pass>
void FreeFermionWavefunction::Amplitude::do_perform_move (const Move &move)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction *>(wf.get());

    if (first_pass) {
        if (wf_->jastrow) {
            m_current_jastrow = wf_->jastrow->compute_jastrow(r);
            if (m_current_jastrow == amplitude_t(0)) {
                m_partial_update_step = get_N_species();
                return;
            }
        } else {
            BOOST_ASSERT(m_current_jastrow == 1);
        }
    }

    for (unsigned int i = (first_pass ? 0 : get_N_species() - m_partial_update_step); i < get_N_species(); ++i) {
        if (m_species_move_in_progress[i]) {
            lw_vector<std::pair<unsigned int, unsigned int>, MAX_MOVE_SIZE> cols;
            for (unsigned int j = 0; j < move.size(); ++j) {
                if (move[j].particle.species == i)
                    cols.push_back(std::make_pair(move[j].particle.index, move[j].destination));
            }
            BOOST_ASSERT(cols.size() != 0);
            m_cmat[i].update_columns(cols, wf_->orbital_def[i]->get_orbitals());
            if (first_pass && m_cmat[i].is_singular()) {
                m_partial_update_step = get_N_species() - i - 1;
                return;
            }
        }
    }
    m_partial_update_step = 0;
}

Big<amplitude_t> FreeFermionWavefunction::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return Big<amplitude_t>();

    Big<amplitude_t> rv(m_current_jastrow);
    for (unsigned int i = 0; i < get_N_species(); ++i)
        rv *= m_cmat[i].get_determinant();
    return rv;
}

void FreeFermionWavefunction::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    for (unsigned int i = 0; i < get_N_species(); ++i) {
        if (m_species_move_in_progress[i])
            m_cmat[i].finish_columns_update();
    }
}

void FreeFermionWavefunction::Amplitude::cancel_move_ (void)
{
    // we use an int here because the condition (i >= 0) would always be true
    // if we used an unsigned int!
    for (int i = get_N_species() - 1 - m_partial_update_step; i >= 0; --i) {
        if (m_species_move_in_progress[i])
            m_cmat[i].cancel_columns_update();
    }
    m_current_jastrow = m_old_jastrow;
    m_partial_update_step = 0;
}

void FreeFermionWavefunction::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    BOOST_ASSERT(species < m_cmat.size());
    m_cmat[species].swap_columns(particle1_index, particle2_index);
}

void FreeFermionWavefunction::Amplitude::reset_ (const PositionArguments &r_)
{
    r = r_;
    m_cmat.resize(0);
    reinitialize();
}

void FreeFermionWavefunction::Amplitude::reinitialize (void)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction *>(wf.get());

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    BOOST_ASSERT(r.get_N_species() == get_N_species());
    BOOST_ASSERT(r.get_N_sites() == wf_->orbital_def[0]->get_N_sites());
    for (unsigned int i = 0; i < get_N_species(); ++i)
        BOOST_ASSERT(r.get_N_filled(i) == wf_->orbital_def[i]->get_N_filled());
#endif

    if (wf_->jastrow) {
        m_current_jastrow = wf_->jastrow->compute_jastrow(r);
    }

    BOOST_ASSERT(m_cmat.size() == 0);
    for (unsigned int j = 0; j < wf_->orbital_def.size(); ++j) {
        const unsigned int N = r.get_N_filled(j);
        Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
        for (unsigned int i = 0; i < N; ++i)
            mat.col(i) = wf_->orbital_def[j]->at_position(r[Particle(i, j)]);
        m_cmat.push_back(CeperleyMatrix<amplitude_t>(mat));
    }
}

boost::shared_ptr<Wavefunction::Amplitude> FreeFermionWavefunction::Amplitude::clone_ (void) const
{
    return boost::make_shared<FreeFermionWavefunction::Amplitude>(*this);
}
