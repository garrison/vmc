#ifndef _SUBSYSTEM_OCCUPATION_NUMBER_PROBABILITY_MEASUREMENT_HPP
#define _SUBSYSTEM_OCCUPATION_NUMBER_PROBABILITY_MEASUREMENT_HPP

#include <vector>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "Measurement.hpp"
#include "StandardWalk.hpp"
#include "Subsystem.hpp"

/**
 * Subsystem occupation number probability measurement
 *
 * Calculates the probability of each possible subsystem occupation number
 * (i.e. how often is a precise, given number of particles of each species
 * found within the subsystem?)
 */
class SubsystemOccupationNumberProbabilityMeasurement : public Measurement<StandardWalk>
{
public:
    SubsystemOccupationNumberProbabilityMeasurement (const boost::shared_ptr<const Subsystem> &subsystem_,
                                                     unsigned int steps_per_measurement)
        : Measurement<StandardWalk>(steps_per_measurement),
          subsystem(subsystem_)
        {
        }

    /**
     * Get the probability so far for a given occupation
     *
     * @param occupation a vector, each item of which represents a subsystem
     * occupation count for the corresponding species
     *
     * @return the probability of the given occupation
     *
     * When iterating through this, pay attention to cache characteristics:
     * it's best to iterate through the first item of the vector, then the
     * second, etc.
     */
    real_t get (const std::vector<unsigned int> &occupation) const
        {
            BOOST_ASSERT(occupation.size() == offsets.size());

            // find the "offset" for the value we are looking for.  this is
            // really just a way of doing dynamic multi-dimensional arrays in
            // C++ without a bunch of overhead.
            unsigned int offset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
                offset += occupation[i] * offsets[i];

            // if the following assertion fails, the occupation must contain an
            // invalid value that is greater than the respective
            // r.get_N_filled(species)
            BOOST_ASSERT(offset < data.size());

            return real_t(data[offset]) / get_measurements_completed();
        }

private:
    void initialize_ (const StandardWalk &walk)
        {
            BOOST_ASSERT(data.size() == 0);

            const PositionArguments &r = walk.get_wavefunction().get_positions();

            // set up the `offsets` and `data` vectors
            offsets.resize(r.get_N_species());
            unsigned int n = 1;
            for (unsigned int i = 0; i < r.get_N_species(); ++i) {
                offsets[i] = n;
                n *= r.get_N_filled(i);
            }
            data.resize(n);
        }

    void measure_ (const StandardWalk &walk)
        {
            const WavefunctionAmplitude &wf = walk.get_wavefunction();

            // calculate the "offset" (see above)
            unsigned int offset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
                offset += do_subsystem_particle_count(wf, i) * offsets[i];

            last_data_ptr = &data[offset];
            ++data[offset];
        }

    void repeat_measurement_ (const StandardWalk &walk)
        {
            (void) walk;
            ++(*last_data_ptr);
        }

    unsigned int do_subsystem_particle_count (const WavefunctionAmplitude &wf, unsigned int species) const
        {
            const PositionArguments &r = wf.get_positions();
            unsigned int rv = 0;
            for (unsigned int i = 0; i < r.get_N_filled(species); ++i) {
                const Particle particle(i, species);
                if (subsystem->position_is_within(r[particle], wf.get_lattice()))
                    ++rv;
            }
            return rv;
        }

    const boost::shared_ptr<const Subsystem> subsystem;
    std::vector<unsigned int> data;
    std::vector<unsigned int> offsets;
    unsigned int *last_data_ptr; // saves us from recalculating the offset on repeat measurements
};

#endif