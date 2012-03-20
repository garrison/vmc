#ifndef _SWAPPED_SYSTEM_HPP
#define _SWAPPED_SYSTEM_HPP

#include <vector>

#include <boost/shared_ptr.hpp>

#include "PositionArguments.hpp"
#include "Subsystem.hpp"
#include "WavefunctionAmplitude.hpp"

class Lattice;

/**
 * Keeps tracked of the swapped wavefunctions
 *
 * phialpha refers to the unswapped wavefunctions, and phibeta refers to the
 * swapped wavefunctions.  This language (and the concept for this class) are
 * from the paper Y. Zhang et. al., PRL 107, 067202 (2011).
 *
 * Any time the phialpha's are updated, this object must be notified
 * immediately after.  At most one particle can be moved in each phialpha at a
 * time.  This is fine, as both RenyiModMeasurement and RenyiSignWalk obey this
 * restriction.
 *
 * @see RenyiModMeasurement
 * @see RenyiSignWalk
 */
class SwappedSystem
{
public:
    SwappedSystem (const boost::shared_ptr<const Subsystem> &subsystem_);

    void initialize (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);

    /**
     * Update the phibetas
     *
     * This must be called immediately after a particle is moved in either or
     * both of the phialpha's.
     *
     * particle1 and particle2 refer to which particle moved in phialpha1 and
     * phialpha2 respectively.  It is possible for one of these to be null,
     * which signifies that no particle was moved in that copy.
     */
    void update (const Particle *particle1, const Particle *particle2,
                 const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);

    /**
     * Finish the update of the phibetas
     */
    void finish_update (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);

    const Subsystem & get_subsystem (void) const
        {
            return *subsystem;
        }

    const WavefunctionAmplitude & get_phibeta1 (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            BOOST_ASSERT(subsystem_particle_counts_match());
            return *phibeta1;
        }

    const WavefunctionAmplitude & get_phibeta2 (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            BOOST_ASSERT(subsystem_particle_counts_match());
            return *phibeta2;
        }

    /**
     *
     */
    bool subsystem_particle_counts_match (void) const;

private:
    void reinitialize_phibetas (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);
    void verify_phibetas (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2) const;
    void swap_positions (PositionArguments &r1, PositionArguments &r2) const;

    boost::shared_ptr<WavefunctionAmplitude> phibeta1, phibeta2; // copy on write
    const boost::shared_ptr<const Subsystem> subsystem;
    std::vector<std::vector<unsigned int> > copy1_subsystem_indices, copy2_subsystem_indices;
    bool phibeta1_dirty, phibeta2_dirty; // helps save time on RenyiSign calculation

    enum NextStep {
        INITIALIZE,
        UPDATE,
        FINISH_UPDATE
    } next_step;
};

/**
 *
 */
extern bool count_subsystem_particle_counts_for_match (const WavefunctionAmplitude &wf1, const WavefunctionAmplitude &wf2,
                                                       const Subsystem &subsystem);

#endif
