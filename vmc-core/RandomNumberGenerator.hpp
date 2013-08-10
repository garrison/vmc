#ifndef _VMC_RANDOM_NUMBER_GENERATOR_HPP
#define _VMC_RANDOM_NUMBER_GENERATOR_HPP

#include <memory>

#include <boost/noncopyable.hpp>

/**
 * Generic interface for a random number generator
 */
class RandomNumberGenerator : boost::noncopyable
{
public:
    typedef unsigned long long rng_seed_t;

    virtual ~RandomNumberGenerator (void)
        {
        }

    /**
     * Returns a random integer in the range [lower_bound, upper_cutoff)
     */
    virtual unsigned int random_small_uint (unsigned int lower_bound, unsigned int upper_cutoff) = 0;

    /**
     * Returns a random integer in the range [0, upper_cutoff)
     */
    unsigned int random_small_uint (unsigned int upper_cutoff)
        {
            return random_small_uint(0, upper_cutoff);
        }

    /**
     * Returns a double in the range [0, 1)
     */
    virtual double random_uniform01 (void) = 0;

    /**
     * Returns a number from the Gaussian distribution centered at 0 with
     * $\sigma=1$
     */
    virtual double random_gaussian (void) = 0;

    /**
     * Returns true if the given name represents a valid random number
     * generator
     */
    static bool name_is_valid (const char *rng_name);

    /**
     * Creates a specific type of RandomNumberGenerator, given by name
     */
    static std::unique_ptr<RandomNumberGenerator> create (const char *rng_name, rng_seed_t seed);
};

#endif
