#include <iostream>

#include "MetropolisSimulation.hpp"
#include "Chain1d.hpp"

int main ()
{
    Chain1d wf(24, 48);
    Chain1dWalk walk(wf);
    MetropolisSimulation<Chain1dWalk> vmc_sim(walk, 16);
    vmc_sim.iterate(16);
    std::cout << "Accepted " << vmc_sim.steps_accepted() << " of " << vmc_sim.steps_completed() << " steps." << std::endl;
}
