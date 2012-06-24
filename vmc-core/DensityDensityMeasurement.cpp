#include "DensityDensityMeasurement.hpp"

real_t DensityDensityMeasurement::get (unsigned int site_index, unsigned int basis_index) const
{
    BOOST_ASSERT(site_index < density_accum.cols());
    BOOST_ASSERT(basis_index < density_accum.rows());
    unsigned int num = density_accum(basis_index, site_index);
    return real_t(num) / denominator - density_squared;
}

void DensityDensityMeasurement::initialize_ (const StandardWalk &walk)
{
    const unsigned int total_sites = walk.get_wavefunction().get_lattice().total_sites();
    BOOST_ASSERT(total_sites > 0);
    const Lattice *lattice = &walk.get_wavefunction().get_lattice();

    const unsigned int basis_indices = lattice->basis_indices;
    density_accum.setZero(basis_indices, total_sites);
    current_step_density_accum.resizeLike(density_accum);

    single_step_denominator = lattice->total_sites() / basis_indices;

    const PositionArguments &r = walk.get_wavefunction().get_positions();

    BOOST_ASSERT(species1 < r.get_N_species());
    BOOST_ASSERT(species2 < r.get_N_species());

    density_squared = real_t(r.get_N_filled(species1) * r.get_N_filled(species2)) / (r.get_N_sites() * r.get_N_sites());
}

void DensityDensityMeasurement::measure_ (const StandardWalk &walk)
{
    const PositionArguments &r = walk.get_wavefunction().get_positions();
    const Lattice *lattice = &walk.get_wavefunction().get_lattice();

    current_step_density_accum.setZero();

    // loop through all pairs of particles
    for (unsigned int i = 0; i < r.get_N_filled(species2); ++i) {
        const LatticeSite site_i(lattice->site_from_index(r[Particle(i, species2)]));
        for (unsigned int j = 0; j < r.get_N_filled(species1); ++j) {
            LatticeSite site_j(lattice->site_from_index(r[Particle(j, species1)]));
            lattice->asm_subtract_site_vector(site_j, site_i.bravais_site());
            lattice->enforce_boundary(site_j);
            ++current_step_density_accum(site_i.basis_index, lattice->site_to_index(site_j));
        }
    }

    repeat_measurement_(walk);
}

void DensityDensityMeasurement::repeat_measurement_ (const StandardWalk &walk)
{
    (void) walk;
    density_accum += current_step_density_accum;
    denominator += single_step_denominator;
}
