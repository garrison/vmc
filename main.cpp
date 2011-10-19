#include <iostream>
#include <vector>
#include <memory>

#include "MetropolisSimulation.hpp"
#include "Chain1d.hpp"

const int N = 200;
const int sample_size = 8;

const int subsystem_upper_bound = 0;

template <typename T>
static inline T square (T v)
{
    return v * v;
}

static std::pair<double, double> stats (const std::vector<double> &v)
// returns mean, stddev of mean
{
    double mean = 0;
    for (unsigned int i = 0; i < v.size(); ++i)
	mean += v[i];
    mean /= v.size();

    double variance = 0;
    for (unsigned int i = 0; i < v.size(); ++i)
	variance += square(v[i] - mean);

    return std::pair<double, double>(mean, sqrt(variance / (v.size() - 1) / v.size()));
}

int main ()
{
#if 1
    rng_class rng(8);
    Chain1d wf(N / 2, N);
    std::vector<boost::shared_ptr<Chain1dContiguousSubsystem> > subsystem;
    std::vector<std::vector<boost::shared_ptr<MetropolisSimulation<Chain1dRenyiWalk, Chain1dRenyiMeasurement> > > > vmc_sim(N);
    for (int i = 0; i <= subsystem_upper_bound; ++i) {
	std::cerr << "Initializing subsystem " << i << std::endl;
	subsystem.push_back(boost::shared_ptr<Chain1dContiguousSubsystem>(new Chain1dContiguousSubsystem(20)));
	for (int j = 0; j < sample_size; ++j) {
	    Chain1dRenyiWalk walk(wf, &*subsystem[i], Chain1dRenyiWalk::SWAPA_MOD, rng);
	    vmc_sim[i].push_back(boost::shared_ptr<MetropolisSimulation<Chain1dRenyiWalk, Chain1dRenyiMeasurement> >(new MetropolisSimulation<Chain1dRenyiWalk, Chain1dRenyiMeasurement>(walk, 12, i * sample_size + j)));
	}
    }
    for (;;) {
	for (int i = 0; i <= subsystem_upper_bound; ++i) {
	    std::cerr << "Iterating on subsystem " << i << std::endl;
	    for (int j = 0; j < sample_size; ++j)
		vmc_sim[i][j]->iterate(12);
	    std::cout << "Current measurement on subsystem (" << i << "): ";
	    for (int j = 0; j < sample_size; ++j)
		std::cout << vmc_sim[i][j]->get_measurement() << " (S" << vmc_sim[i][j]->get_walk().get_N_subsystem() << " " << (static_cast<double>(vmc_sim[i][j]->steps_accepted()) / vmc_sim[i][j]->steps_completed() * 100) << "%)  ";
	    std::cout << std::endl;

	    // calculate and display stats
	    std::vector<double> v;
	    for (int j = 0; j < sample_size; ++j)
		v.push_back(vmc_sim[i][j]->get_measurement());
	    std::pair<double, double> s = stats(v);
	    std::cerr << s.first << " \\pm " << s.second << std::endl;
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
