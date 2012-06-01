/**
 * Main file for RVB ( = fully projected BCS wave function of spinons at half-filling ) stuff.
 */

#include <iostream>
#include <vector>
#include <list>
#include <memory>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "DensityDensityMeasurement.hpp"
#include "FilledOrbitals.hpp"
#include "SimpleSubsystem.hpp"
#include "HypercubicLattice.hpp"
#include "RVBWavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-filling.hpp"
#include "allowed-momentum.hpp"
#include "lowest-momenta.hpp"

const unsigned int seed = 13;

const unsigned int DIMENSION = 2;
const unsigned int lattice_length = 6;

int main ()
{
    // initialize random number generator
    rng_class rng(seed);
    
    // set up a lattice
    lw_vector<int, MAX_DIMENSION> lattice_dimensions(DIMENSION);
    for (unsigned int i = 0; i < DIMENSION; ++i)
        lattice_dimensions[i] = lattice_length;
    boost::shared_ptr<const HypercubicLattice<DIMENSION> > lattice(new HypercubicLattice<DIMENSION>(lattice_dimensions));
    const unsigned int M = lattice->total_sites()/2;
        
    // set up initial particle positions at random
    std::vector<std::vector<unsigned int> > vv(2);
    vv[0] = some_random_filling(M, *lattice, rng);
    vv[1] = some_random_filling(M, *lattice, rng);
    PositionArguments r(vv, lattice->total_sites());
    
    // RVB wave function with the "pair wave function" in the determinant a Gaussian
    std::vector<complex_t> phi;
    for (unsigned int i = 0; i < lattice->total_sites(); ++i) phi.push_back(exp(-1.0*i*i));
    boost::shared_ptr<WavefunctionAmplitude> wf(new RVBWavefunctionAmplitude(r, lattice, phi));
    
    // try different initial particle positions until a non-zero amplitude is found
    const bool success = search_for_filling_with_nonzero_amplitude(*wf, rng);
    if (!success) {
        std::cerr << "could not find filling with non-zero amplitude in a reasonable amount of time" << std::endl;
        return 1;
    }
    
    // set up three monte carlo simulations and perform some number of steps to
    // get each to equilibrium
    StandardWalk walk(wf);
    boost::shared_ptr<DensityDensityMeasurement> density_measurement(new DensityDensityMeasurement(10, 0, 0));
    MetropolisSimulation<StandardWalk> sim(walk, density_measurement, 5000, rng());
        
    // continue iterating on each simulation, outputting results periodically
    for (unsigned int i = 0; i < 20; ++i) {
        sim.iterate(5000);
        std::cerr << "density-density " << (100.0 * sim.steps_accepted() / sim.steps_completed()) << "%\t";
        for (unsigned int i = 0; i < lattice->total_sites(); ++i)
            std::cerr << "  " << density_measurement->get(i);
        std::cerr << std::endl;        
    }
}
