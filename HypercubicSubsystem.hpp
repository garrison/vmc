#ifndef _HYPERCUBIC_SUBSYSTEM_HPP
#define _HYPERCUBIC_SUBSYSTEM_HPP

#include <boost/array.hpp>

#include "Subsystem.hpp"
#include "HypercubicLattice.hpp"

template <std::size_t DIM>
class HypercubicSubsystem : public Subsystem
{
public:
    HypercubicSubsystem (unsigned int subsystem_length_)
	{
	    for (unsigned int i = 0; i < DIM; ++i)
		subsystem_length[i] = subsystem_length_;
	}

    HypercubicSubsystem (const boost::array<unsigned int, DIM> &subsystem_length_)
	: subsystem_length(subsystem_length_)
	{
	}

    bool particle_is_within (unsigned int site_index, const Lattice &lattice_) const
	{
	    BOOST_ASSERT(lattice_makes_sense(lattice_));
	    const HypercubicLattice<DIM> *lattice = dynamic_cast<const HypercubicLattice<DIM> *>(&lattice_);
	    BOOST_ASSERT(lattice != 0);

	    typename HypercubicLattice<DIM>::Site site(lattice->site_from_index(site_index));
	    for (unsigned int i = 0; i < DIM; ++i) {
		BOOST_ASSERT(site[i] >= 0);
		if (site[i] >= (int) subsystem_length[i])
		    return false;
	    }
	    return true;
	}

    bool lattice_makes_sense (const Lattice &lattice) const
	{
	    return bool(dynamic_cast<const HypercubicLattice<DIM> *>(&lattice));
	}

private:
    boost::array<unsigned int, DIM> subsystem_length;

    // disable default constructor
    HypercubicSubsystem (void);
};

#endif
