#include <boost/assert.hpp>
#include <boost/make_shared.hpp>

#include "vmc-math-utils.hpp"
#include "DBLWavefunctionAmplitude.hpp"

DBLWavefunctionAmplitude::DBLWavefunctionAmplitude (const PositionArguments &r_, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_1, const boost::shared_ptr<const OrbitalDefinitions> &orbital_def_2, real_t d1_exponent_, real_t d2_exponent_)
    : WavefunctionAmplitude(r_, orbital_def_1->get_lattice_ptr()),
      orbital_def1(orbital_def_1),
      orbital_def2(orbital_def_2),
      d1_exponent(d1_exponent_),
      d2_exponent(d2_exponent_),
      m_partial_update_step(0)
{
    BOOST_ASSERT(r.get_N_species() == 1);
    BOOST_ASSERT(r.get_N_sites() == orbital_def1->get_N_sites());
    BOOST_ASSERT(r.get_N_sites() == orbital_def2->get_N_sites());
    BOOST_ASSERT(r.get_N_filled(0) == orbital_def1->get_N_filled());
    BOOST_ASSERT(r.get_N_filled(0) == orbital_def2->get_N_filled());
    BOOST_ASSERT(orbital_def1->get_lattice_ptr() == orbital_def2->get_lattice_ptr());

    reinitialize();
}

void DBLWavefunctionAmplitude::perform_move_ (const Move &move)
{
    // we require that m_partial_update_step == 0 between moves; otherwise, psi_() will
    // return zero when it shouldn't.
    BOOST_ASSERT(m_partial_update_step == 0);

    m_current_move = move;

    do_perform_move<true>(move);
}

// this function exists solely to get around a bug in which clang++ fails to
// compile if we perform this operation directly
static inline void set_column_from_orbitals (Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> &mat, unsigned int column, const OrbitalDefinitions &orbitals, unsigned int destination)
{
    mat.col(column) = orbitals.at_position(destination);
}

// During perform_move(), we want to short-circuit out as soon as either of the
// determinants is zero, as we already know what psi_() will be.  However,
// doing so means we must later make a second pass to finish the (possibly)
// remaining determinant update if finish_update() is called.  This templated
// function allows us to use the same code in both passes, with only minor
// differences.
template <bool first_pass>
void DBLWavefunctionAmplitude::do_perform_move (const Move &move)
{
    if (!first_pass && m_partial_update_step == 0)
        return;

    const unsigned int N = orbital_def1->get_N_filled();

    lw_vector<unsigned int, MAX_MOVE_SIZE> c;
    for (unsigned int i = 0; i < move.size(); ++i)
        c.push_back(move[i].particle.index);

    switch (first_pass ? 2 : m_partial_update_step)
    {
    case 2:
        {
            Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> d2_cols(N, move.size());
            for (unsigned int i = 0; i < move.size(); ++i)
                set_column_from_orbitals(d2_cols, i, *orbital_def2, move[i].destination);
            cmat2.update_columns(c, d2_cols);
        }
        if (first_pass && cmat2.get_determinant() == amplitude_t(0)) {
            m_partial_update_step = 1;
            return;
        }

    case 1:
        {
            Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> d1_cols(N, move.size());
            for (unsigned int i = 0; i < move.size(); ++i)
                set_column_from_orbitals(d1_cols, i, *orbital_def1, move[i].destination);
            cmat1.update_columns(c, d1_cols);
        }
    }

    if (!first_pass)
        m_partial_update_step = 0;
}

amplitude_t DBLWavefunctionAmplitude::psi_ (void) const
{
    if (m_partial_update_step != 0)
        return amplitude_t(0);

    // fixme: we could cache or precalculate this ... but i doubt it would make
    // much difference really
    return (complex_pow(cmat1.get_determinant(), d1_exponent)
            * complex_pow(cmat2.get_determinant(), d2_exponent));
}

void DBLWavefunctionAmplitude::finish_move_ (void)
{
    do_perform_move<false>(m_current_move);

    cmat1.finish_columns_update();
    cmat2.finish_columns_update();
}

void DBLWavefunctionAmplitude::cancel_move_ (void)
{
    switch (m_partial_update_step) {
    case 0:
        cmat1.cancel_columns_update();
    case 1:
        cmat2.cancel_columns_update();
    }

    m_partial_update_step = 0;
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
