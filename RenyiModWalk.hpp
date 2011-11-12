#ifndef _RENYI_MOD_WALK_HPP
#define _RENYI_MOD_WALK_HPP

#include <vector>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>

#include "vmc-typedefs.hpp"
#include "Subsystem.hpp"
#include "SwappedSystem.hpp"
#include "WavefunctionAmplitude.hpp"

class RenyiModWalk
{
public:
    RenyiModWalk (const boost::shared_ptr<WavefunctionAmplitude> &wf, const std::vector<boost::shared_ptr<const Subsystem> > &subsystems, rng_class &rng);
    probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
    void accept_transition (void);

    const WavefunctionAmplitude & get_phialpha1 (void) const
	{
	    return *phialpha1;
	}

    const WavefunctionAmplitude & get_phialpha2 (void) const
	{
	    return *phialpha2;
	}

    const WavefunctionAmplitude & get_phibeta1 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index]->get_phibeta1();
	}

    const WavefunctionAmplitude & get_phibeta2 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index]->get_phibeta2();
	}

    unsigned int get_N_subsystem1 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index]->get_N_subsystem1();
	}

    unsigned int get_N_subsystem2 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index]->get_N_subsystem2();
	}

    unsigned int get_subsystem_array_size (void) const
	{
	    return swapped_system.size();
	}

private:
    boost::shared_ptr<WavefunctionAmplitude> phialpha1, phialpha2;
    std::vector<boost::shared_ptr<SwappedSystem> > swapped_system;
    unsigned int transition_copy_in_progress;
    unsigned int chosen_particle;

    // disable the default constructor
    RenyiModWalk (void);
};

#endif
