#ifndef _VMC_BASE_SWAP_POSSIBLE_WALK_HPP
#define _VMC_BASE_SWAP_POSSIBLE_WALK_HPP

#include <memory>
#include <utility>

#include <boost/assert.hpp>

#include "vmc-typedefs.hpp"
#include "Subsystem.hpp"
#include "SwappedSystem.hpp"
#include "Wavefunction.hpp"
#include "Walk.hpp"

class RandomNumberGenerator;

/**
 * Base class for Renyi walks that require the same number/types of particles
 * in the subsystem for each copy (i.e., a swap is always possible)
 *
 * See Y. Zhang et. al., PRL 107, 067202 (2011) for explanation
 *
 * @see RenyiSignWalk
 * @see RenyiModPossibleWalk
 */
class BaseSwapPossibleWalk : public Walk<probability_t>
{
public:
    /**
     * Constructor.
     *
     * It is essential that both wfa1 and wfa2 have the same number of
     * particles in the subsystem.
     */
    BaseSwapPossibleWalk (std::unique_ptr<Wavefunction<amplitude_t>::Amplitude> wfa1, std::unique_ptr<Wavefunction<amplitude_t>::Amplitude> wfa2, const std::shared_ptr<const Subsystem> &subsystem, bool update_swapped_system_before_accepting_=true);

    const Wavefunction<amplitude_t>::Amplitude & get_phialpha1 (void) const
        {
            return *phialpha1;
        }

    const Wavefunction<amplitude_t>::Amplitude & get_phialpha2 (void) const
        {
            return *phialpha2;
        }

    const Wavefunction<amplitude_t>::Amplitude & get_phibeta1 (void) const
        {
            return swapped_system.get_phibeta1();
        }

    const Wavefunction<amplitude_t>::Amplitude & get_phibeta2 (void) const
        {
            return swapped_system.get_phibeta2();
        }

private:
    virtual probability_t compute_probability_ratio_of_random_transition (RandomNumberGenerator &rng) override final;
    virtual void accept_transition (void) override;
    virtual void reject_transition (void) override;
    virtual void check_for_numerical_error (void) const override;

    virtual probability_t probability_ratio (amplitude_t phialpha1_ratio, amplitude_t phialpha2_ratio, amplitude_t phibeta1_ratio, amplitude_t phibeta2_ratio) const = 0;

    static unsigned int count_subsystem_sites (const Subsystem &subsystem, const Lattice &lattice);

    std::unique_ptr<Wavefunction<amplitude_t>::Amplitude> phialpha1, phialpha2;
    SwappedSystem swapped_system;
    Particle chosen_particle_A, chosen_particle_B;
    const Particle *chosen_particle1, *chosen_particle2;
    unsigned int N_subsystem_sites; // remains constant after initialization
    bool update_swapped_system_before_accepting; // remains constant after initialization
    bool autoreject_in_progress;
    bool transition_in_progress;
};

#endif
