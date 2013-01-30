#include <boost/noncopyable.hpp>

#include "OperatorMeasurement.hpp"

/**
 * Manages a temporary move on a wave function.
 *
 * When it goes out of scope, the move is automatically cancelled.
 *
 * NEVER allow an update while this is in scope, or create a second one on a
 * wavefunction when it is in scope.  Keep it around for as little time as
 * possible.
 */
class TemporaryMove : boost::noncopyable
{
public:
    TemporaryMove (const Wavefunction::Amplitude &wfa_, const Move &move)
        : wfa(wfa_)
        {
            const_cast<Wavefunction::Amplitude &>(wfa).perform_move(move);
        }

    ~TemporaryMove (void)
        {
            const_cast<Wavefunction::Amplitude &>(wfa).cancel_move();
        }

    const Wavefunction::Amplitude &wfa;
};

void OperatorMeasurement::initialize_ (const StandardWalk &walk)
{
    (void) walk;
    BOOST_ASSERT(&walk.get_wavefunctionamplitude().get_lattice() == m_operator.lattice.get());
}

void OperatorMeasurement::measure_ (const StandardWalk &walk)
{
    const Wavefunction::Amplitude &wfa = walk.get_wavefunctionamplitude();
    const PositionArguments &r = wfa.get_positions();
    const Lattice &lattice = wfa.get_lattice();

    amplitude_t meas = 0;

    const amplitude_t old_psi = wfa.psi();

    // we only iterate if doing a sum, and even then we only want to iterate
    // over BraivaisSite's
    const unsigned int n_iterations = sum ? lattice.total_sites() / lattice.basis_indices : 1;

    // FIXME: need a faster way of iterating the lattice
    for (unsigned int i = 0; i < n_iterations; ++i) {
        const LatticeSite site_offset(lattice.site_from_index(i));
        // we only want to iterate over BravaisSite's
        BOOST_ASSERT(site_offset.basis_index == 0);

        phase_t phase = 1;

        Move move;
        for (unsigned int j = 0; j < m_operator.hopv.size(); ++j) {
            phase_t srcphase;
            const unsigned int species = m_operator.hopv[j].get_species();
            LatticeSite src(m_operator.hopv[j].get_source());
            lattice.asm_add_site_vector(src, site_offset.bravais_site());
            srcphase = lattice.enforce_boundary(src, bcs);
            if (srcphase == phase_t(0))
                goto current_measurement_is_zero;
            const int particle_index = r.particle_index_at_position(lattice.site_to_index(src), species);
            if (particle_index < 0)
                goto current_measurement_is_zero;
            if (m_operator.hopv[j].get_source() != m_operator.hopv[j].get_destination()) {
                LatticeSite dest(m_operator.hopv[j].get_destination());
                lattice.asm_add_site_vector(dest, site_offset.bravais_site());
                phase *= lattice.enforce_boundary(dest, bcs) / srcphase;
                if (phase == phase_t(0))
                    goto current_measurement_is_zero;
                const unsigned int dest_index = lattice.site_to_index(dest);
                if (r.is_occupied(dest_index, species))
                    goto current_measurement_is_zero;
                move.push_back(SingleParticleMove(Particle(particle_index, species), dest_index));
            }
        }

        // now perform the move (if necessary)
        if (move.size() != 0) {
            TemporaryMove temp_move(wfa, move);
            // fixme: check logic of multiplying by phase (c.f. above), as
            // well as logic of source and destination
            meas += std::conj(wfa.psi() * phase / old_psi);
        } else {
            meas += 1;
        }

    current_measurement_is_zero: ;
    }

    most_recent_value = meas;
    estimate.add_value(most_recent_value);
}

void OperatorMeasurement::repeat_measurement_ (const StandardWalk &walk)
{
    (void) walk;
    estimate.add_value(most_recent_value);
}
