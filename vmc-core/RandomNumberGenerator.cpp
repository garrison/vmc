#include <cstring>

#include <boost/assert.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>

#include "RandomNumberGenerator.hpp"

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
    unsigned int random_small_uint (unsigned int lower_bound, unsigned int upper_cutoff)
        {
            BOOST_ASSERT(lower_bound < upper_cutoff);
            boost::uniform_smallint<> distribution(lower_bound, upper_cutoff - 1);
            boost::variate_generator<T&, boost::uniform_smallint<> > generator(rng, distribution);
            return generator();
        }

    double random_uniform01 (void)
        {
            return uniform01_distribution();
        }

private:
    T rng;

    // see <http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html> for why we
    // must declare this here
    boost::uniform_01<T> uniform01_distribution;
};

std::auto_ptr<RandomNumberGenerator> RandomNumberGenerator::create (const char *rng_name, RandomNumberGenerator::rng_seed_t seed)
{
    if (std::strcmp(rng_name, "boost::mt19937") == 0)
        return std::auto_ptr<RandomNumberGenerator>(new RandomNumberGeneratorBoostImpl<boost::mt19937>(seed));

    BOOST_ASSERT(false); // invalid rng specified
    return std::auto_ptr<RandomNumberGenerator>(); // suppress non-void return warning
}
