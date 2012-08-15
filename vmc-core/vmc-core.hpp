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
    HighlevelSimulation (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice);

    void iterate (unsigned int sweeps)
        {
            sim->iterate(sweeps);
        }

    std::string output (void) const;

private:
    boost::shared_ptr<RandomNumberGenerator> rng;
    boost::shared_ptr<BaseMetropolisSimulation> sim;
    std::list<boost::shared_ptr<BaseMeasurement> > measurements;
    std::string walk_type;
};

#endif
