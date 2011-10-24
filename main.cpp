#include <iostream>
#include <vector>
#include <memory>

#include "MetropolisSimulation.hpp"
#include "Chain1d.hpp"

const int N = 200;
const int sample_size = 3;

const int subsystem_upper_bound = 0;

const int seed = 43; // 19

typedef Chain1dRenyiModMeasurement CurrentMeasurement;
typedef Chain1dRenyiModWalk CurrentWalk;

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
    std::vector<boost::shared_ptr<Chain1dContiguousSubsystem> > subsystem;
    std::vector<std::vector<boost::shared_ptr<MetropolisSimulation<CurrentWalk, CurrentMeasurement> > > > vmc_sim(N);
    for (int i = 0; i <= subsystem_upper_bound; ++i) {
	std::cerr << "Initializing subsystem " << i << std::endl;
	subsystem.push_back(boost::shared_ptr<Chain1dContiguousSubsystem>(new Chain1dContiguousSubsystem(100)));
	for (int j = 0; j < sample_size; ++j) {
	    CurrentWalk walk(wf, &*subsystem[i], rng);
	    vmc_sim[i].push_back(boost::shared_ptr<MetropolisSimulation<CurrentWalk, CurrentMeasurement> >(new MetropolisSimulation<CurrentWalk, CurrentMeasurement>(walk, 12, i * sample_size + j)));
	}
    }
    for (;;) {
	for (int i = 0; i <= subsystem_upper_bound; ++i) {
	    std::cerr << "Iterating on subsystem " << i << std::endl;
	    for (int j = 0; j < sample_size; ++j)
		vmc_sim[i][j]->iterate(12);
	    std::cout << "Current measurement on subsystem (" << i << "): ";
	    for (int j = 0; j < sample_size; ++j)
		std::cout << vmc_sim[i][j]->get_measurement() << " (S" << vmc_sim[i][j]->get_walk().get_N_subsystem1() << "," << vmc_sim[i][j]->get_walk().get_N_subsystem2() << " " << (static_cast<double>(vmc_sim[i][j]->steps_accepted()) / vmc_sim[i][j]->steps_completed() * 100) << "%)  ";
	    std::cout << std::endl;

	    // calculate and display stats
	    std::vector<CurrentMeasurement::measurement_value_t> v;
	    for (int j = 0; j < sample_size; ++j)
		v.push_back(vmc_sim[i][j]->get_measurement());
	    std::pair<CurrentMeasurement::measurement_value_t, CurrentMeasurement::measurement_value_t> s = stats(v);
	    std::cerr << s.first << " \\pm " << s.second << std::endl;
	}
    }
}
