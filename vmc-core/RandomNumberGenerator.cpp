#include <cstring>

#include <boost/assert.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>

#include "RandomNumberGenerator.hpp"
#include "vmcstd.hpp"

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
            BOOST_ASSERT(lower_bound < upper_cutoff);
            if (lower_bound == upper_cutoff - 1)
                return lower_bound; // there's only one possibility, so we don't need randomness
            boost::uniform_smallint<> distribution(lower_bound, upper_cutoff - 1);
            boost::variate_generator<T&, boost::uniform_smallint<> > generator(rng, distribution);
            return generator();
        }

    virtual double random_uniform01 (void) override
        {
            return uniform01_distribution();
        }

private:
    T rng;

    // see <http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html> for why we
    // must declare this here
    boost::uniform_01<T> uniform01_distribution;
};

static const char * const rng_names[] = {
    "boost::mt19937",
    "boost::lagged_fibonacci607",
    0 // a null pointer marks the end of the array
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
    BOOST_ASSERT(name_is_valid(rng_name));
    if (std::strcmp(rng_name, "boost::mt19937") == 0)
        return vmcstd::make_unique<RandomNumberGeneratorBoostImpl<boost::mt19937> >(seed);
    if (std::strcmp(rng_name, "boost::lagged_fibonacci607") == 0)
        return vmcstd::make_unique<RandomNumberGeneratorBoostImpl<boost::lagged_fibonacci607> >(seed);

    BOOST_ASSERT(false); // invalid rng specified
    return std::unique_ptr<RandomNumberGenerator>(); // suppress non-void return warning
}
