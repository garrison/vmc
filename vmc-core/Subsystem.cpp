#include "Subsystem.hpp"
#include "WavefunctionAmplitude.hpp"

unsigned int count_N_subsystem (const WavefunctionAmplitude &wf, const Subsystem &subsystem)
{
    BOOST_ASSERT(subsystem.lattice_makes_sense(wf.get_lattice()));
    const PositionArguments &r = wf.get_positions();
    unsigned int rv = 0;
    for (unsigned int i = 0; i < r.size(); ++i) {
        if (subsystem.position_is_within(r[i], wf.get_lattice()))
            ++rv;
    }
    return rv;
}
