#include <boost/noncopyable.hpp>

#include "GreenMeasurement.hpp"

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
    TemporaryMove (const WavefunctionAmplitude &wf_, const Move &move)
        : wf(wf_)
        {
            const_cast<WavefunctionAmplitude &>(wf).perform_move(move);
        }

    ~TemporaryMove (void)
        {
            const_cast<WavefunctionAmplitude &>(wf).cancel_move();
        }

    const WavefunctionAmplitude &wf;
};

void GreenMeasurement::initialize_ (const StandardWalk &walk)
{
    const unsigned int total_sites = walk.get_wavefunction().get_lattice().total_sites();
    BOOST_ASSERT(total_sites > 0);
    const Lattice *lattice = &walk.get_wavefunction().get_lattice();

    const unsigned int basis_indices = lattice->basis_indices;
    green_accum.setZero(basis_indices, total_sites);
    current_step_green_accum.resizeLike(green_accum);

    single_step_denominator = lattice->total_sites() / basis_indices;
}

void GreenMeasurement::measure_ (const StandardWalk &walk)
{
    const WavefunctionAmplitude &wf = walk.get_wavefunction();
    const PositionArguments &r = wf.get_positions();
    const Lattice *lattice = &wf.get_lattice();

    current_step_green_accum.setZero();

    const amplitude_t old_psi = wf.psi();

    // loop through all (particle, empty site) pairs
    for (unsigned int i = 0; i < r.get_N_filled(species); ++i) {
        const Particle particle(i, species);
        const LatticeSite site_i(lattice->site_from_index(r[particle]));

        // amplitude of starting and ending on same site
        current_step_green_accum(site_i.basis_index, 0) += amplitude_t(1);

        // loop through all empty sites
        for (unsigned int j = 0; j < r.get_N_sites(); ++j) {
            if (r.is_occupied(j, species))
                continue;

            Move move;
            move.push_back(SingleParticleMove(particle, j));
            TemporaryMove temp_move(wf, move);

            LatticeSite site_j(lattice->site_from_index(j));
            lattice->asm_subtract_site_vector(site_j, site_i.bravais_site());
            const phase_t phase = lattice->enforce_boundary(site_j);
            // fixme: check logic of multiplying by phase
            current_step_green_accum(site_i.basis_index, lattice->site_to_index(site_j)) += std::conj(wf.psi() * phase / old_psi);
        }
    }

    repeat_measurement_(walk);
}

void GreenMeasurement::repeat_measurement_ (const StandardWalk &walk)
{
    (void) walk;
    green_accum += current_step_green_accum;
    denominator += single_step_denominator;
}
