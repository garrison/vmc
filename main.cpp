#include <iostream>

#include "MetropolisSimulation.hpp"
#include "Chain1d.hpp"

int main ()
{
    Chain1d wf(4, 8);
    MetropolisSimulation<Chain1d> vmc_sim(wf, 12);
    vmc_sim.iterate(8);
    std::cout << "Accepted " << vmc_sim.steps_accepted() << " of " << vmc_sim.steps_completed() << " steps." << std::endl;
}
