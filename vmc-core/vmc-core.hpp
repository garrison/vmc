#ifndef _VMC_CORE_HPP
#define _VMC_CORE_HPP

#include <list>

#include <boost/shared_ptr.hpp>

#include "MetropolisSimulation.hpp"
#include "Measurement.hpp"
#include "Lattice.hpp"

MetropolisSimulation * create_simulation (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice, const std::list<boost::shared_ptr<BaseMeasurement> > &measurements, unsigned int equilibrium_steps);

#endif
