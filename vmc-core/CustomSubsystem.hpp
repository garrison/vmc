#ifndef _VMC_CUSTOM_SUBSYSTEM_HPP
#define _VMC_CUSTOM_SUBSYSTEM_HPP

#include <boost/dynamic_bitset.hpp>

#include "Subsystem.hpp"

/**
 * Subsystem containing some custom set of lattice sites
 */
class CustomSubsystem : public Subsystem
{
public:
    /**
     * @param site_status_ a bitset of the size of the lattice, denoting each site's presence in the subsystem
     */
    explicit CustomSubsystem (const boost::dynamic_bitset<> &site_status_)
        : site_status(site_status_)
        {
        }

    virtual bool position_is_within (unsigned int site_index, const Lattice &lattice) const override;

    virtual bool lattice_makes_sense (const Lattice &lattice) const override;

    const boost::dynamic_bitset<> site_status;
};

#endif
