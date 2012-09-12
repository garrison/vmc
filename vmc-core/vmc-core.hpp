#ifndef _VMC_CORE_HPP
#define _VMC_CORE_HPP

#include <list>
#include <memory>

#include <boost/shared_ptr.hpp>

#include "MetropolisSimulation.hpp"
#include "Measurement.hpp"
#include "Lattice.hpp"

std::auto_ptr<Walk> create_walk_from_json (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice, std::auto_ptr<RandomNumberGenerator> &rng);

#endif
