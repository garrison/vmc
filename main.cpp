#include <iostream>
#include <vector>
#include <memory>

#include "MetropolisSimulation.hpp"
#include "Chain1d.hpp"

const int N = 200;
const int sample_size = 3;

int main ()
{
#if 1
    rng_class rng(3);
    Chain1d wf(N / 2, N);
    std::vector<boost::shared_ptr<Chain1dContiguousSubsystem> > subsystem;
    std::vector<std::vector<boost::shared_ptr<MetropolisSimulation<Chain1dRenyiWalk> > > > vmc_sim(N);
    for (int i = 0; i < 4; ++i) {
	std::cerr << "Initializing subsystem " << i << std::endl;
	subsystem.push_back(boost::shared_ptr<Chain1dContiguousSubsystem>(new Chain1dContiguousSubsystem(i)));
	Chain1dRenyiWalk walk(wf, &*subsystem[i], Chain1dRenyiWalk::SWAPA_MOD, rng);
	for (int j = 0; j < sample_size; ++j)
	    vmc_sim[i].push_back(boost::shared_ptr<MetropolisSimulation<Chain1dRenyiWalk> >(new MetropolisSimulation<Chain1dRenyiWalk>(walk, 12, i * sample_size + j)));
    }
    for (;;) {
	for (int i = 0; i < N; ++i) {
	    std::cerr << "Iterating on subsystem " << i << std::endl;
	    for (int j = 0; j < sample_size; ++j)
		vmc_sim[i][j]->iterate(12);
	    std::cout << "Current measurement on subsystem " << i << ": " << std::endl;
	}
    }
#endif

#if 0
    std::cout << "Accepted " << vmc_sim.steps_accepted() << " of " << vmc_sim.steps_completed() << " steps." << std::endl;
    vmc_sim.get_walk();
#endif

#if 0
    static rng_class generator(4);
    boost::uniform_smallint<> integer_distribution(0, 10);
    boost::variate_generator<rng_class&, boost::uniform_smallint<> > gen(generator, integer_distribution);

    unsigned int N = 20;
    std::vector<int> v(N);
    for (int i = 0; i < 200000; ++i) {
	for (unsigned int j = 0; j < N; ++j)
	    v[j] = gen();
	Chain1dContiguousSubsystem subsystem(gen());
	std::pair<int, int> p = fermion_partition(v, subsystem);
	int particles = p.first, parity = p.second;
	for (std::vector<int>::const_iterator i = v.begin(); i != v.end(); ++i)
	    std::cout << *i << ' ';
	std::cout << "\nparticles: " << particles << "\nparity: " << parity << std::endl;
    }
#endif
}
