/**
 * This file is designed to give a glimpse of the C++ API, without being as
 * verbose as vmc-core.cpp.
 */

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
#include "DBLWavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-filling.hpp"
#include "allowed-momentum.hpp"
#include "lowest-momenta.hpp"

const unsigned int seed = 13;

const unsigned int DIMENSION = 1;
const unsigned int lattice_length = 40;
const unsigned int F = 20;

int main ()
{
    // initialize random number generator
    rng_class rng(seed);

    // set up a lattice
    boost::array<int, DIMENSION> lattice_dimensions;
    for (unsigned int i = 0; i < DIMENSION; ++i)
        lattice_dimensions[i] = lattice_length;
    boost::shared_ptr<const HypercubicLattice<DIMENSION> > lattice(new HypercubicLattice<DIMENSION>(lattice_dimensions));

    // set up initial particle positions at random
    std::vector<std::vector<unsigned int> > vv;
    vv.push_back(some_random_filling<DIMENSION>(F, *lattice, rng));
    PositionArguments r(vv, lattice->total_sites());

    // set up the boundary conditions and orbitals.  this will give a DBL with
    // identical fermi seas.
    HypercubicLattice<DIMENSION>::BoundaryConditions boundary_conditions;
    for (unsigned int i = 0; i < DIMENSION; ++i)
        boundary_conditions[i] = periodic_bc;
    boost::shared_ptr<const OrbitalDefinitions> orbitals(new FilledOrbitals<DIMENSION>(lowest_momenta(*lattice, boundary_conditions, F), lattice, boundary_conditions));
    boost::shared_ptr<WavefunctionAmplitude> wf(new DBLWavefunctionAmplitude(r, orbitals, orbitals, 1.0, 1.0));

    // try different initial particle positions until a non-zero amplitude is found
    const bool success = search_for_filling_with_nonzero_amplitude<DIMENSION>(*wf, *lattice, rng);
    if (!success) {
        std::cerr << "could not find filling with non-zero amplitude in a reasonable amount of time" << std::endl;
        return 1;
    }

    // set up three monte carlo simulations and perform some number of steps to
    // get each to equilibrium
    StandardWalk walk(wf);
    boost::shared_ptr<DensityDensityMeasurement<DIMENSION> > density_measurement(new DensityDensityMeasurement<DIMENSION>(1, 0, 0));
    MetropolisSimulation<StandardWalk> sim(walk, density_measurement, 5000, rng());

    std::list<boost::shared_ptr<Measurement<RenyiModWalk> > > mod_measurements;
    mod_measurements.push_back(boost::make_shared<RenyiModMeasurement>(boost::make_shared<SimpleSubsystem<DIMENSION> >(4), 50));

    RenyiModWalk mod_walk(wf, wf);
    MetropolisSimulation<RenyiModWalk> mod_sim(mod_walk, mod_measurements, 5000, rng());

    boost::shared_ptr<Subsystem> subsystem(new SimpleSubsystem<DIMENSION>(4));
    RenyiSignWalk sign_walk(wf, wf, subsystem);
    boost::shared_ptr<RenyiSignMeasurement> sign_measurement(new RenyiSignMeasurement);
    MetropolisSimulation<RenyiSignWalk> sign_sim(sign_walk, sign_measurement, 5000, rng());

    // continue iterating on each simulation, outputting results periodically
    for (unsigned int i = 0; i < 20; ++i) {
        sim.iterate(5000);
        std::cerr << "density-density " << (100.0 * sim.steps_accepted() / sim.steps_completed()) << "%\t";
        for (unsigned int i = 0; i < lattice->total_sites(); ++i)
            std::cerr << "  " << density_measurement->get(i);
        std::cerr << std::endl;

        mod_sim.iterate(5000);
        std::cerr << "swap,mod " << (100.0 * mod_sim.steps_accepted() / mod_sim.steps_completed()) << "%\t" << double(dynamic_cast<RenyiModMeasurement *>(&**mod_measurements.begin())->get()) << std::endl;

        sign_sim.iterate(5000);
        std::cerr << "swap,sign " << (100.0 * sign_sim.steps_accepted() / sign_sim.steps_completed()) << "%\t" << sign_measurement->get() << std::endl;
    }
}
