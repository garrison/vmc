#ifndef _VMC_SUBSYSTEM_OCCUPATION_NUMBER_PROBABILITY_MEASUREMENT_HPP
#define _VMC_SUBSYSTEM_OCCUPATION_NUMBER_PROBABILITY_MEASUREMENT_HPP

#include <vector>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "Measurement.hpp"
#include "BlockedEstimate.hpp"
#include "StandardWalk.hpp"
#include "Wavefunction.hpp"
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
    SubsystemOccupationNumberProbabilityMeasurement (unsigned int steps_per_measurement, const boost::shared_ptr<const Subsystem> &subsystem_)
        : Measurement<StandardWalk>(steps_per_measurement),
          subsystem(subsystem_)
        {
        }

    /**
     * Get the probability estimate so far for a given occupation
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
    const BlockedEstimate<unsigned int> & get_estimate (const std::vector<unsigned int> &occupation) const
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
            BOOST_ASSERT(offset < estimate.size());

            return estimate[offset];
        }

    const std::vector<unsigned int> & get_bounds (void) const
        {
            return bounds;
        }

private:
    void initialize_ (const StandardWalk &walk)
        {
            BOOST_ASSERT(estimate.size() == 0);

            const PositionArguments &r = walk.get_wavefunctionamplitude().get_positions();

            // set up the `offsets`, `bounds`, and `estimate` vectors
            offsets.resize(r.get_N_species());
            bounds.resize(r.get_N_species());
            unsigned int n = 1;
            for (unsigned int i = 0; i < r.get_N_species(); ++i) {
                offsets[i] = n;
                // we need (N_filled + 1) slots since the number in the
                // subsystem will be in the range 0 .. N_filled
                n *= r.get_N_filled(i) + 1;
                bounds[i] = r.get_N_filled(i);
            }
            estimate.resize(n);
        }

    void measure_ (const StandardWalk &walk)
        {
            const Wavefunction::Amplitude &wfa = walk.get_wavefunctionamplitude();

            // calculate the "offset" (see above)
            unsigned int offset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
                offset += do_subsystem_particle_count(wfa, i) * offsets[i];
            BOOST_ASSERT(offset < estimate.size());

            last_offset = offset;
            repeat_measurement_(walk);
        }

    void repeat_measurement_ (const StandardWalk &walk)
        {
            (void) walk;
            // tally the value 1 for the current occupation, and 0 for each
            // other possible occupation
            estimate[last_offset].add_value(1);
            for (unsigned int i = 0; i < estimate.size(); ++i) {
                if (i != last_offset)
                    estimate[i].add_value(0);
            }
        }

    unsigned int do_subsystem_particle_count (const Wavefunction::Amplitude &wfa, unsigned int species) const
        {
            const PositionArguments &r = wfa.get_positions();
            unsigned int rv = 0;
            for (unsigned int i = 0; i < r.get_N_filled(species); ++i) {
                const Particle particle(i, species);
                if (subsystem->position_is_within(r[particle], wfa.get_lattice()))
                    ++rv;
            }
            return rv;
        }

    const boost::shared_ptr<const Subsystem> subsystem;
    std::vector<BlockedEstimate<unsigned int> > estimate;
    std::vector<unsigned int> offsets, bounds;
    unsigned int last_offset; // saves us from recalculating the offset on repeat measurements
};

#endif
