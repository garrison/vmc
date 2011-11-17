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
#include "HypercubicSubsystem.hpp"
#include "FreeFermionWavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-combination.hpp"

const unsigned int seed = 13;

const unsigned int N = 40;
const unsigned int F = N / 2;

int main ()
{
    rng_class rng(seed);

    std::vector<unsigned int> v;
    random_combination(v, F, N, rng);

    PositionArguments r(v, N);
    boost::shared_ptr<WavefunctionAmplitude> wf(new FreeFermionWavefunctionAmplitude(r));

    StandardWalk walk(wf);
    boost::shared_ptr<DensityDensityMeasurement<HypercubicLattice<1> > > density_measurement(new DensityDensityMeasurement<HypercubicLattice<1> >);
    MetropolisSimulation<StandardWalk> sim(walk, density_measurement, 8, rng());

    std::list<boost::shared_ptr<Measurement<RenyiModWalk> > > mod_measurements;
    mod_measurements.push_back(boost::make_shared<RenyiModMeasurement>(boost::make_shared<HypercubicSubsystem<1> >(4)));

    RenyiModWalk mod_walk(wf, rng);
    MetropolisSimulation<RenyiModWalk> mod_sim(mod_walk, mod_measurements, 8, rng());

    boost::shared_ptr<Subsystem> subsystem(new HypercubicSubsystem<1>(4));
    RenyiSignWalk sign_walk(wf, subsystem, rng);
    boost::shared_ptr<RenyiSignMeasurement> sign_measurement(new RenyiSignMeasurement);
    MetropolisSimulation<RenyiSignWalk> sign_sim(sign_walk, sign_measurement, 8, rng());

    for (unsigned int i = 0; i < 20; ++i) {
	sim.iterate(12);
	std::cerr << "density-density " << (100.0 * sim.steps_accepted() / sim.steps_completed()) << "%\t";
	for (unsigned int i = 0; i < N; ++i)
	    std::cerr << "  " << density_measurement->get(i);
	std::cerr << std::endl;

	mod_sim.iterate(12);
	std::cerr << "swap,mod " << (100.0 * mod_sim.steps_accepted() / mod_sim.steps_completed()) << "%\t" << double(dynamic_cast<RenyiModMeasurement *>(&**mod_measurements.begin())->get()) << std::endl;

	sign_sim.iterate(12);
	std::cerr << "swap,sign " << (100.0 * sign_sim.steps_accepted() / sign_sim.steps_completed()) << "%\t" << sign_measurement->get() << std::endl;
    }
}
