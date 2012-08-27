#ifndef _VMC_CORE_HPP
#define _VMC_CORE_HPP

#include <string>
#include <list>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "MetropolisSimulation.hpp"
#include "Measurement.hpp"
#include "Lattice.hpp"

class HighlevelSimulation : boost::noncopyable
{
public:
    HighlevelSimulation (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice, const std::list<boost::shared_ptr<BaseMeasurement> > &measurements);

    void iterate (unsigned int sweeps)
        {
            sim->iterate(sweeps);
        }

    std::string output (void) const;

    unsigned int steps_completed (void) const
        {
            return sim->steps_completed();
        }

    unsigned int steps_accepted (void) const
        {
            return sim->steps_accepted();
        }

    unsigned int steps_fully_rejected (void) const
        {
            return sim->steps_fully_rejected();
        }

private:
    boost::shared_ptr<MetropolisSimulation> sim;
    std::string walk_type;
};

#endif
