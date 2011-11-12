#include <iostream>
#include <vector>
#include <memory>

#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "RenyiModMeasurement.hpp"
#include "RenyiModWalk.hpp"
#include "RenyiSignMeasurement.hpp"
#include "RenyiSignWalk.hpp"
#include "NullMeasurement.hpp"
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
    MetropolisSimulation<StandardWalk, NullMeasurement<StandardWalk> > sim(walk, 8, rng());
    sim.iterate(12);

    boost::shared_ptr<const Subsystem> subsystem(new HypercubicSubsystem<1>(4));
    std::vector<boost::shared_ptr<const Subsystem> > subsystems;
    subsystems.push_back(subsystem);

    RenyiModWalk mod_walk(wf, subsystems, rng);
    MetropolisSimulation<RenyiModWalk, RenyiModMeasurement> mod_sim(mod_walk, 8, rng());
    mod_sim.iterate(12);

    RenyiSignWalk sign_walk(wf, subsystem, rng);
    MetropolisSimulation<RenyiSignWalk, RenyiSignMeasurement> sign_sim(sign_walk, 8, rng());
    sign_sim.iterate(12);
}
