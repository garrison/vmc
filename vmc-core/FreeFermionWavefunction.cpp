#include <boost/assert.hpp>
#include <boost/cast.hpp>

#include "FreeFermionWavefunction.hpp"

template <typename AmplitudeType>
FreeFermionWavefunction<AmplitudeType>::FreeFermionWavefunction (const std::vector<std::shared_ptr<const OrbitalDefinitions<AmplitudeType> > > &orbital_def_, const std::shared_ptr<const JastrowFactor<AmplitudeType> > &jastrow_)
    : Wavefunction<AmplitudeType>(orbital_def_[0]->get_lattice_ptr()),
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

template <typename AmplitudeType>
FreeFermionWavefunction<AmplitudeType>::Amplitude::Amplitude (const std::shared_ptr<const FreeFermionWavefunction> &wf_, const PositionArguments &r_)
    : Wavefunction<AmplitudeType>::Amplitude(wf_, r_),
      m_current_jastrow(1),
      m_partial_update_step(0),
      m_species_move_in_progress(wf_->get_N_species())
{
    reinitialize();
}

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::perform_move_ (const Move &move)
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
template <typename AmplitudeType>
template <bool first_pass>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::do_perform_move (const Move &move)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction<AmplitudeType> *>(this->wf.get());

    if (first_pass) {
        if (wf_->jastrow) {
            m_current_jastrow = wf_->jastrow->compute_jastrow(this->r);
            if (m_current_jastrow.is_zero()) {
                m_partial_update_step = get_N_species();
                return;
            }
        } else {
            BOOST_ASSERT(m_current_jastrow.get_base() == 1. && m_current_jastrow.get_exponent() == 0.);
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

template <typename AmplitudeType>
Big<AmplitudeType> FreeFermionWavefunction<AmplitudeType>::Amplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return Big<AmplitudeType>();

    Big<AmplitudeType> rv(m_current_jastrow);
    for (unsigned int i = 0; i < get_N_species(); ++i)
        rv *= m_cmat[i].get_determinant();
    return rv;
}

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    for (unsigned int i = 0; i < get_N_species(); ++i) {
        if (m_species_move_in_progress[i])
            m_cmat[i].finish_columns_update();
    }
}

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::cancel_move_ (void)
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

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::swap_particles_ (unsigned int particle1_index, unsigned int particle2_index, unsigned int species)
{
    BOOST_ASSERT(species < m_cmat.size());
    m_cmat[species].swap_columns(particle1_index, particle2_index);

    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction<AmplitudeType> *>(this->wf.get());

    if (wf_->jastrow) {
        m_current_jastrow = wf_->jastrow->compute_jastrow(this->r);
    }
}

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::reset_ (const PositionArguments &r_)
{
    this->r = r_;
    m_cmat.resize(0);
    reinitialize();
}

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::reinitialize (void)
{
    const FreeFermionWavefunction *wf_ = boost::polymorphic_downcast<const FreeFermionWavefunction *>(this->wf.get());

#if !defined(BOOST_DISABLE_ASSERTS) && !defined(NDEBUG)
    BOOST_ASSERT(this->r.get_N_species() == get_N_species());
    BOOST_ASSERT(this->r.get_N_sites() == wf_->orbital_def[0]->get_N_sites());
    for (unsigned int i = 0; i < get_N_species(); ++i)
        BOOST_ASSERT(this->r.get_N_filled(i) == wf_->orbital_def[i]->get_N_filled());
#endif

    if (wf_->jastrow) {
        m_current_jastrow = wf_->jastrow->compute_jastrow(this->r);
    }

    BOOST_ASSERT(m_cmat.size() == 0);
    for (unsigned int j = 0; j < wf_->orbital_def.size(); ++j) {
        const unsigned int N = this->r.get_N_filled(j);
        Eigen::Matrix<AmplitudeType, Eigen::Dynamic, Eigen::Dynamic> mat(N, N);
        for (unsigned int i = 0; i < N; ++i)
            mat.col(i) = wf_->orbital_def[j]->at_position(this->r[Particle(i, j)]);
        m_cmat.push_back(CeperleyMatrix<AmplitudeType>(mat));
    }
}

template <typename AmplitudeType>
void FreeFermionWavefunction<AmplitudeType>::Amplitude::check_for_numerical_error (void) const
{
    for (unsigned int i = 0; i < get_N_species(); ++i)
        m_cmat[i].check_for_numerical_error();
}

template <typename AmplitudeType>
std::shared_ptr<typename Wavefunction<AmplitudeType>::Amplitude> FreeFermionWavefunction<AmplitudeType>::Amplitude::clone_ (void) const
{
    return std::make_shared<FreeFermionWavefunction<AmplitudeType>::Amplitude>(*this);
}

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) template class FreeFermionWavefunction<type>
#include "vmc-supported-types.hpp"
