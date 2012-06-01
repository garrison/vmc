#ifndef _SIMPLE_SUBSYSTEM_HPP
#define _SIMPLE_SUBSYSTEM_HPP

#include <boost/array.hpp>
#include <boost/assert.hpp>

#include "lw_vector.hpp"
#include "Subsystem.hpp"
#include "Lattice.hpp"

/**
 * Represents any subsystem that is a parallelpiped aligned with the lattice's
 * primitive vectors
 */
template <std::size_t DIM>
class SimpleSubsystem : public Subsystem
{
public:
    explicit SimpleSubsystem (const lw_vector<unsigned int, MAX_DIMENSION> &subsystem_length_)
        : subsystem_length(subsystem_length_)
        {
        }

    bool position_is_within (unsigned int site_index, const Lattice &lattice) const
        {
            BOOST_ASSERT(lattice_makes_sense(lattice));

            const LatticeSite site(lattice.site_from_index(site_index));
            for (unsigned int i = 0; i < lattice.n_dimensions(); ++i) {
                BOOST_ASSERT(site[i] >= 0);
                if (site[i] >= (int) subsystem_length[i])
                    return false;
            }
            return true;
        }

    bool lattice_makes_sense (const Lattice &lattice) const
        {
            if (lattice.n_dimensions() != subsystem_length.size())
                return false;
            for (unsigned int i = 0; i < lattice.n_dimensions(); ++i) {
                if (lattice.dimensions[i] < (int) subsystem_length[i])
                    return false;
            }
            return true;
        }

private:
    const lw_vector<unsigned int, MAX_DIMENSION> subsystem_length;
};

#endif
