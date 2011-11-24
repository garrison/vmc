#include <iostream>
#include <vector>
#include <list>
#include <memory>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "RenyiModMeasurement.hpp"
#include "RenyiModWalk.hpp"
#include "RenyiSignMeasurement.hpp"
#include "RenyiSignWalk.hpp"
#include "DensityDensityMeasurement.hpp"
#include "FilledOrbitals.hpp"
#include "SimpleSubsystem.hpp"
#include "HypercubicLattice.hpp"
#include "FreeFermionWavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-combination.hpp"
#include "allowed-momentum.hpp"
#include "lowest-momenta.hpp"
#include "array-util.hpp"

const unsigned int seed = 13;

const unsigned int DIMENSION = 1;
const unsigned int lattice_length = 40;
const unsigned int F = 20;

int main ()
{
    rng_class rng(seed);

    boost::shared_ptr<const HypercubicLattice<DIMENSION> > lattice(new HypercubicLattice<DIMENSION>(make_array<int, DIMENSION>(lattice_length)));

    std::vector<unsigned int> v;
    random_combination(v, F, lattice->total_sites(), rng);
    PositionArguments r(v, lattice->total_sites());

    HypercubicLattice<DIMENSION>::BoundaryConditions boundary_conditions = make_array<BoundaryCondition>(periodic_bc);
    boost::shared_ptr<const OrbitalDefinitions> orbitals(new FilledOrbitals<DIMENSION>(lowest_momenta(*lattice, boundary_conditions, F), lattice, boundary_conditions));
    boost::shared_ptr<WavefunctionAmplitude> wf(new FreeFermionWavefunctionAmplitude(r, orbitals));

    StandardWalk walk(wf);
    boost::shared_ptr<DensityDensityMeasurement<DIMENSION> > density_measurement(new DensityDensityMeasurement<DIMENSION>);
    MetropolisSimulation<StandardWalk> sim(walk, density_measurement, 8, rng());

    std::list<boost::shared_ptr<Measurement<RenyiModWalk> > > mod_measurements;
    mod_measurements.push_back(boost::make_shared<RenyiModMeasurement>(boost::make_shared<SimpleSubsystem<DIMENSION> >(4)));

    RenyiModWalk mod_walk(wf, rng);
    MetropolisSimulation<RenyiModWalk> mod_sim(mod_walk, mod_measurements, 8, rng());

    boost::shared_ptr<Subsystem> subsystem(new SimpleSubsystem<DIMENSION>(4));
    RenyiSignWalk sign_walk(wf, subsystem, rng);
    boost::shared_ptr<RenyiSignMeasurement> sign_measurement(new RenyiSignMeasurement);
    MetropolisSimulation<RenyiSignWalk> sign_sim(sign_walk, sign_measurement, 8, rng());

    for (unsigned int i = 0; i < 20; ++i) {
        sim.iterate(12);
        std::cerr << "density-density " << (100.0 * sim.steps_accepted() / sim.steps_completed()) << "%\t";
        for (unsigned int i = 0; i < lattice->total_sites(); ++i)
            std::cerr << "  " << density_measurement->get(i);
        std::cerr << std::endl;

        mod_sim.iterate(12);
        std::cerr << "swap,mod " << (100.0 * mod_sim.steps_accepted() / mod_sim.steps_completed()) << "%\t" << double(dynamic_cast<RenyiModMeasurement *>(&**mod_measurements.begin())->get()) << std::endl;

        sign_sim.iterate(12);
        std::cerr << "swap,sign " << (100.0 * sign_sim.steps_accepted() / sign_sim.steps_completed()) << "%\t" << sign_measurement->get() << std::endl;
    }
}
