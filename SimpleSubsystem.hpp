#ifndef _SIMPLE_SUBSYSTEM_HPP
#define _SIMPLE_SUBSYSTEM_HPP

#include <boost/array.hpp>
#include <boost/assert.hpp>

#include "Subsystem.hpp"
#include "NDLattice.hpp"

/**
 * Represents any subsystem that is a parallelpiped aligned with the lattice's
 * primitive vectors
 */
template <std::size_t DIM>
class SimpleSubsystem : public Subsystem
{
public:
    SimpleSubsystem (unsigned int subsystem_length_)
        {
            for (unsigned int i = 0; i < DIM; ++i)
                subsystem_length[i] = subsystem_length_;
        }

    SimpleSubsystem (const boost::array<unsigned int, DIM> &subsystem_length_)
        : subsystem_length(subsystem_length_)
        {
        }

    bool position_is_within (unsigned int site_index, const Lattice &lattice_) const
        {
            BOOST_ASSERT(lattice_makes_sense(lattice_));
            const NDLattice<DIM> *lattice = dynamic_cast<const NDLattice<DIM> *>(&lattice_);
            BOOST_ASSERT(lattice != 0);

            typename NDLattice<DIM>::Site site(lattice->site_from_index(site_index));
            for (unsigned int i = 0; i < DIM; ++i) {
                BOOST_ASSERT(site[i] >= 0);
                if (site[i] >= (int) subsystem_length[i])
                    return false;
            }
            return true;
        }

    bool lattice_makes_sense (const Lattice &lattice) const
        {
            const NDLattice<DIM> *ndlattice = dynamic_cast<const NDLattice<DIM> *>(&lattice);
            if (!ndlattice)
                return false;
            for (unsigned int i = 0; i < DIM; ++i) {
                if (ndlattice->length[i] < (int) subsystem_length[i])
                    return false;
            }
            return true;
        }

private:
    boost::array<unsigned int, DIM> subsystem_length;
};

#endif
