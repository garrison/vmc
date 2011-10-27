#include <iostream>
#include <vector>
#include <memory>

#include "MetropolisSimulation.hpp"
#include "Chain1d.hpp"

const unsigned int N = 10;
const unsigned int sample_size = 8;

const int seed = 56;

template <typename T>
static inline T square (T v)
{
    return v * v;
}

template <typename T>
static std::pair<T, T> stats (const std::vector<T> &v)
// returns mean, stddev of mean
{
    T mean = 0;
    for (unsigned int i = 0; i < v.size(); ++i)
	mean += v[i];
    mean /= v.size();

    T variance = 0;
    for (unsigned int i = 0; i < v.size(); ++i)
	variance += square(v[i] - mean);

    double denominator = (v.size() - 1) * v.size();
    return std::pair<T, T>(mean, sqrt(variance / denominator));
}

int main ()
{
    rng_class rng(seed);
    Chain1d wf(N / 2, N);
    std::vector<boost::shared_ptr<Chain1dContiguousSubsystem> > subsystems;
    std::vector<boost::shared_ptr<MetropolisSimulation<Chain1dRenyiModWalk, Chain1dRenyiModMeasurement> > > vmc_sim;
    for (unsigned int i = 0; i <= N / 2 + 1; ++i) {
	subsystems.push_back(boost::shared_ptr<Chain1dContiguousSubsystem>(new Chain1dContiguousSubsystem(i)));
    }
    for (unsigned int j = 0; j < sample_size; ++j) {
	Chain1dRenyiModWalk walk(wf, &*subsystems[0], rng);
	vmc_sim.push_back(boost::shared_ptr<MetropolisSimulation<Chain1dRenyiModWalk, Chain1dRenyiModMeasurement> >(new MetropolisSimulation<Chain1dRenyiModWalk, Chain1dRenyiModMeasurement>(walk, 12, rng())));
    }
    for (;;) {
	for (unsigned int j = 0; j < sample_size; ++j)
	    vmc_sim[j]->iterate(12);
	for (unsigned int j = 0; j < sample_size; ++j)
	    std::cout << vmc_sim[j]->get_measurement() << " (S" << vmc_sim[j]->get_walk().get_N_subsystem1() << "," << vmc_sim[j]->get_walk().get_N_subsystem2() << " " << (static_cast<double>(vmc_sim[j]->steps_accepted()) / vmc_sim[j]->steps_completed() * 100) << "%)  ";
	std::cout << std::endl;

	// calculate and display stats
	std::vector<Chain1dRenyiModMeasurement::measurement_value_t> v;
	for (unsigned int j = 0; j < sample_size; ++j)
	    v.push_back(vmc_sim[j]->get_measurement());
	std::pair<Chain1dRenyiModMeasurement::measurement_value_t, Chain1dRenyiModMeasurement::measurement_value_t> s = stats(v);
//	std::cout << "N" << N << " S" << subsystem_length << "  ";
	std::cout << s.first << " \\pm " << s.second << std::endl;
    }
}
