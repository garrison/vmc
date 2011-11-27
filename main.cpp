#include <iostream>
#include <vector>
#include <list>
#include <memory>
#include <exception>
#include <cstring>

#include <json/json.h>
#include <boost/array.hpp>
#include <boost/assert.hpp>
#include <boost/cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "RenyiModMeasurement.hpp"
#include "RenyiModWalk.hpp"
#include "RenyiSignMeasurement.hpp"
#include "RenyiSignWalk.hpp"
#include "DensityDensityMeasurement.hpp"
#include "FilledOrbitals.hpp"
#include "SimpleSubsystem.hpp"
#include "HypercubicLattice.hpp"
#include "FreeFermionWavefunctionAmplitude.hpp"
#include "PositionArguments.hpp"
#include "random-combination.hpp"
#include "allowed-momentum.hpp"
#include "lowest-momenta.hpp"

// fixme: HypercubicLattice references should instead refer to NDLattice once
// json specifies the momenta (and we no longer need lowest_momenta()

class ParseError : public std::exception
{
public:
    ParseError (void)
        : error_message(default_error_message)
        {
        }

    ParseError (const char *error_message_)
        : error_message(error_message_)
        {
        }

    virtual const char * what (void) const throw()
        {
            return error_message;
        }

private:
    static const char *default_error_message;
    const char *error_message;
};

const char *ParseError::default_error_message = "json input error";

static inline void ensure_object (const Json::Value &jsonvalue)
{
    if (!jsonvalue.isObject())
        throw ParseError("object expected");
}

static inline void ensure_array (const Json::Value &jsonvalue)
{
    if (!jsonvalue.isArray())
        throw ParseError("array expected");
}

static void ensure_required (const Json::Value &jsonvalue, const char * const keys[])
{
    BOOST_ASSERT(jsonvalue.isObject());
    for (const char * const *key = keys; *key != NULL; ++key) {
        if (!jsonvalue.isMember(*key))
            throw ParseError("required keys not all given");
    }
}

static void ensure_only (const Json::Value &jsonvalue, const char * const keys[])
{
    BOOST_ASSERT(jsonvalue.isObject());
    for (Json::Value::const_iterator i = jsonvalue.begin(), e = jsonvalue.end(); i != e; ++i) {
        const char * const *key = keys;
        for (; *key != NULL; ++key) {
            if (strcmp(i.memberName(), *key) == 0)
                break;
        }
        if (key == NULL)
            throw ParseError("too many keys provided");
    }
}

template <unsigned int DIM>
static int do_simulation (const Json::Value &json_input, rng_class &rng);

int main ()
{
    // take json input and perform a simulation

    Json::Value json_input;
    std::cin >> json_input;

    try {
        ensure_object(json_input);
        const char * const json_input_required[] = { "rng", "system", NULL };
        ensure_required(json_input, json_input_required);
        ensure_only(json_input, json_input_required);

        // initialize random number generator
        unsigned int seed;
        const Json::Value &json_rng = json_input["rng"];
        ensure_object(json_rng);
        const char * const json_rng_allowed[] = { "seed", NULL };
        ensure_only(json_rng, json_rng_allowed);
        if (json_rng.isMember("seed")) {
            if (!json_rng["seed"].isIntegral())
                throw ParseError("seed must be correct data type");
            seed = json_rng["seed"].asUInt();
        } else {
            throw ParseError("seed must be given");
        }
        rng_class rng(seed);

        // begin setting up the physical system
        const Json::Value &json_system = json_input["system"];
        ensure_object(json_system);
        const char * const json_system_required[] = { "lattice", NULL };
        ensure_required(json_system, json_system_required);
        ensure_only(json_system, json_system_required);

        // begin setting up the lattice
        const Json::Value &json_lattice = json_system["lattice"];
        ensure_object(json_lattice);
        const char * const json_lattice_required[] = { "size", NULL };
        ensure_required(json_lattice, json_lattice_required);
        ensure_only(json_lattice, json_lattice_required);

        // determine the lattice size/dimension
        const Json::Value &json_lattice_size = json_lattice["size"];
        ensure_array(json_lattice_size);
        unsigned int ndimensions = json_lattice_size.size();
        for (unsigned int i = 0; i < ndimensions; ++i) {
            if (!(json_lattice_size[i].isIntegral() && json_lattice_size[i].asInt() > 0))
                throw ParseError("lattice dimensions must be positive integers");
        }

        // dispatch the remainder of the simulation based on the number of
        // dimensions in the system
        switch (ndimensions) {
        case 1:
            return do_simulation<1>(json_input, rng);
        case 2:
            return do_simulation<2>(json_input, rng);
        default:
            throw ParseError("lattice given has a number of dimensions that is not supported by this build");
        }
    } catch (ParseError e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}

template <unsigned int DIM>
static int do_simulation (const Json::Value &json_input, rng_class &rng)
{
    const Json::Value &json_lattice_size = json_input["system"]["lattice"]["size"];
    boost::array<int, DIM> lattice_size_array;
    for (unsigned int i = 0; i < DIM; ++i)
        lattice_size_array[i] = json_lattice_size[i].asInt();
    boost::shared_ptr<const HypercubicLattice<DIM> > lattice(new HypercubicLattice<DIM>(lattice_size_array));

    // example API usage

    const unsigned int F = 20;

    std::vector<unsigned int> v;
    random_combination(v, F, lattice->total_sites(), rng);
    PositionArguments r(v, lattice->total_sites());

    typename HypercubicLattice<DIM>::BoundaryConditions boundary_conditions;
    for (unsigned int i = 0; i < DIM; ++i)
        boundary_conditions[i] = periodic_bc;

    boost::shared_ptr<const OrbitalDefinitions> orbitals(new FilledOrbitals<DIM>(lowest_momenta(*lattice, boundary_conditions, F), lattice, boundary_conditions));
    boost::shared_ptr<WavefunctionAmplitude> wf(new FreeFermionWavefunctionAmplitude(r, orbitals));

    StandardWalk walk(wf);
    boost::shared_ptr<DensityDensityMeasurement<DIM> > density_measurement(new DensityDensityMeasurement<DIM>);
    MetropolisSimulation<StandardWalk> sim(walk, density_measurement, 8, rng());

    std::list<boost::shared_ptr<Measurement<RenyiModWalk> > > mod_measurements;
    mod_measurements.push_back(boost::make_shared<RenyiModMeasurement>(boost::make_shared<SimpleSubsystem<DIM> >(4)));

    RenyiModWalk mod_walk(wf, rng);
    MetropolisSimulation<RenyiModWalk> mod_sim(mod_walk, mod_measurements, 8, rng());

    boost::shared_ptr<Subsystem> subsystem(new SimpleSubsystem<DIM>(4));
    RenyiSignWalk sign_walk(wf, subsystem, rng);
    boost::shared_ptr<RenyiSignMeasurement> sign_measurement(new RenyiSignMeasurement);
    MetropolisSimulation<RenyiSignWalk> sign_sim(sign_walk, sign_measurement, 8, rng());

    for (unsigned int i = 0; i < 20; ++i) {
        sim.iterate(12);
        std::cerr << "density-density " << (100.0 * sim.steps_accepted() / sim.steps_completed()) << "%\t";
        for (unsigned int i = 0; i < lattice->total_sites(); ++i)
            std::cerr << "  " << density_measurement->get(i);
        std::cerr << std::endl;

        mod_sim.iterate(12);
        std::cerr << "swap,mod " << (100.0 * mod_sim.steps_accepted() / mod_sim.steps_completed()) << "%\t" << double(boost::polymorphic_downcast<RenyiModMeasurement *>(&**mod_measurements.begin())->get()) << std::endl;

        sign_sim.iterate(12);
        std::cerr << "swap,sign " << (100.0 * sign_sim.steps_accepted() / sign_sim.steps_completed()) << "%\t" << sign_measurement->get() << std::endl;
    }

    return 0;
}
