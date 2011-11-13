#ifndef _SWAPPED_SYSTEM_HPP
#define _SWAPPED_SYSTEM_HPP

#include <vector>

#include <boost/shared_ptr.hpp>

#include "PositionArguments.hpp"
#include "Subsystem.hpp"
#include "WavefunctionAmplitude.hpp"

class Lattice;

class SwappedSystem
// this class must be notified any time the phialpha's are changed
{
public:
    SwappedSystem (const boost::shared_ptr<const Subsystem> &subsystem_);

    void initialize (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);

    void update (int index1, int index2, const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);

    void finish_update (void);

    const Subsystem & get_subsystem (void) const
	{
	    return *subsystem;
	}

    const WavefunctionAmplitude & get_phibeta1 (void) const
	{
	    BOOST_ASSERT(next_step != INITIALIZE);
	    BOOST_ASSERT(get_N_subsystem1() == get_N_subsystem2());
	    return *phibeta1;
	}

    const WavefunctionAmplitude & get_phibeta2 (void) const
	{
	    BOOST_ASSERT(next_step != INITIALIZE);
	    BOOST_ASSERT(get_N_subsystem1() == get_N_subsystem2());
	    return *phibeta2;
	}

    unsigned int get_N_subsystem1 (void) const
	{
	    return copy1_subsystem_indices.size();
	}

    unsigned int get_N_subsystem2 (void) const
	{
	    return copy2_subsystem_indices.size();
	}

private:
    void reinitialize_phibetas (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2);
    void verify_phibetas (const WavefunctionAmplitude &phialpha1, const WavefunctionAmplitude &phialpha2) const;
    void swap_positions (PositionArguments &r1, PositionArguments &r2) const;

    // disable default constructor
    SwappedSystem (void);

    boost::shared_ptr<WavefunctionAmplitude> phibeta1, phibeta2; // copy on write
    const boost::shared_ptr<const Subsystem> subsystem;
    std::vector<unsigned int> copy1_subsystem_indices, copy2_subsystem_indices;
    bool phibeta1_dirty, phibeta2_dirty; // helps save time on RenyiSign calculation

    enum NextStep {
	INITIALIZE,
	UPDATE,
	FINISH_UPDATE
    } next_step;
};

#endif
