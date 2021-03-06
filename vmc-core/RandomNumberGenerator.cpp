#include <cassert>
#include <cstring>
#include <random>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>

#include "RandomNumberGenerator.hpp"
#include "vmcstd.hpp"

template <class RandomEngine>
class RandomNumberGeneratorStdImpl : public RandomNumberGenerator
{
public:
    RandomNumberGeneratorStdImpl (RandomNumberGenerator::rng_seed_t seed)
        : engine(seed)
        {
        }

protected:
    virtual unsigned int random_small_uint (unsigned int lower_bound, unsigned int upper_cutoff) override
        {
            assert(lower_bound < upper_cutoff);
            if (lower_bound == upper_cutoff - 1)
                return lower_bound; // there's only one possibility, so we don't need randomness
            std::uniform_int_distribution<unsigned int> distribution(lower_bound, upper_cutoff - 1);
            return distribution(engine);
        }

    virtual int random_sign (void) override
        {
            return std::bernoulli_distribution()(engine) ? 1 : -1;
        }

    virtual double random_uniform01 (void) override
        {
            return std::uniform_real_distribution<>()(engine);
        }

    virtual double random_gaussian (void) override
        {
            return std::normal_distribution<>()(engine);
        }

private:
    RandomEngine engine;
};

template <class T>
class RandomNumberGeneratorBoostImpl : public RandomNumberGenerator
{
public:
    RandomNumberGeneratorBoostImpl (RandomNumberGenerator::rng_seed_t seed)
        : rng(seed),
          uniform01_distribution(rng)
        {
        }

protected:
    virtual unsigned int random_small_uint (unsigned int lower_bound, unsigned int upper_cutoff) override
        {
            assert(lower_bound < upper_cutoff);
            if (lower_bound == upper_cutoff - 1)
                return lower_bound; // there's only one possibility, so we don't need randomness
            boost::uniform_smallint<> distribution(lower_bound, upper_cutoff - 1);
            boost::variate_generator<T&, boost::uniform_smallint<> > generator(rng, distribution);
            return generator();
        }

    virtual int random_sign (void) override
        {
            return (random_small_uint(0, 2) << 1) - 1;
        }

    virtual double random_uniform01 (void) override
        {
            return uniform01_distribution();
        }

    virtual double random_gaussian (void) override
        {
            boost::normal_distribution<> distribution;
            boost::variate_generator<T&, boost::normal_distribution<> > generator(rng, distribution);
            return generator();
        }

private:
    T rng;

    // see <http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html> for why we
    // must declare this here
    boost::uniform_01<T> uniform01_distribution;
};

static const char * const rng_names[] = {
    "std::mt19937",
    "boost::mt19937",
    "boost::lagged_fibonacci607",
    nullptr // marks the end of the array
};

bool RandomNumberGenerator::name_is_valid (const char *rng_name)
{
    for (const char * const * current_rng_name = rng_names; *current_rng_name; ++current_rng_name) {
        if (std::strcmp(rng_name, *current_rng_name) == 0)
            return true;
    }
    return false;
}

std::unique_ptr<RandomNumberGenerator> RandomNumberGenerator::create (const char *rng_name, RandomNumberGenerator::rng_seed_t seed)
{
    assert(name_is_valid(rng_name));
    if (std::strcmp(rng_name, "std::mt19937") == 0)
        return vmcstd::make_unique<RandomNumberGeneratorStdImpl<std::mt19937> >(seed);
    if (std::strcmp(rng_name, "boost::mt19937") == 0)
        return vmcstd::make_unique<RandomNumberGeneratorBoostImpl<boost::mt19937> >(seed);
    if (std::strcmp(rng_name, "boost::lagged_fibonacci607") == 0)
        return vmcstd::make_unique<RandomNumberGeneratorBoostImpl<boost::lagged_fibonacci607> >(seed);

    assert(false); // invalid rng specified
    return std::unique_ptr<RandomNumberGenerator>(); // suppress non-void return warning
}
