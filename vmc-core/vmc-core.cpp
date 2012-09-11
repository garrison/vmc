/*
 * vmc-core.cpp: everything necessary to set up a simulation based on json
 * input, run the simulation, and return its results in json format.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <list>
#include <exception>
#include <cstring>
#include <sstream>
#include <string>

#include <json/json.h>
#include <boost/assert.hpp>
#include <boost/cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <Eigen/Dense>

#include "RandomNumberGenerator.hpp"
#include "Lattice.hpp"
#include "PositionArguments.hpp"
#include "OrbitalDefinitions.hpp"
#include "FreeFermionWavefunction.hpp"
#include "DBLWavefunction.hpp"
#include "DMetalWavefunction.hpp"
#include "RVBWavefunction.hpp"
#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "OperatorMeasurement.hpp"
#include "SubsystemOccupationNumberProbabilityMeasurement.hpp"
#include "RenyiModPossibleWalk.hpp"
#include "RenyiModPossibleMeasurement.hpp"
#include "RenyiSignWalk.hpp"
#include "RenyiSignMeasurement.hpp"
#include "Subsystem.hpp"
#include "SimpleSubsystem.hpp"
#include "vmc-core.hpp"

class ParseError : public std::exception
{
public:
    ParseError (void)
        : error_message(default_error_message)
        {
        }

    // we only store a pointer to the error message, so it must be something
    // that won't be destructed as the exception propagates
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

static void ensure_array (const Json::Value &jsonvalue, unsigned int array_length)
{
    if (!jsonvalue.isArray())
        throw ParseError("array expected");
    if (jsonvalue.size() != array_length)
        throw ParseError("array is not the correct size");
}

static inline void ensure_string (const Json::Value &jsonvalue)
{
    if (!jsonvalue.isString())
        throw ParseError("string expected");
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
            if (std::strcmp(i.memberName(), *key) == 0)
                break;
        }
        if (*key == NULL)
            throw ParseError("too many keys provided");
    }
}

static void ensure_object_with_type_field_as_string (const Json::Value &jsonvalue)
{
    ensure_object(jsonvalue);
    const char * const json_type_required[] = { "type", NULL };
    ensure_required(jsonvalue, json_type_required);
    ensure_string(jsonvalue["type"]);
}

static complex_t parse_complex (const Json::Value &complex)
{
    ensure_array(complex, 2);
    if (!complex[0].isNumeric() || !complex[1].isNumeric())
        throw ParseError("expecting a complex number represented as a json array of two numbers");

    return complex_t(complex[0].asDouble(), complex[1].asDouble());
}

static double json_get_double (const Json::Value &jsonvalue, const char *key, double default_value)
{
    if (!jsonvalue.isMember(key))
        return default_value;

    const Json::Value & jsondouble = jsonvalue[key];
    if (!jsondouble.isNumeric())
        throw ParseError("number expected");
    return jsondouble.asDouble();
}

static boost::shared_ptr<const OrbitalDefinitions> parse_json_orbitals_from_definitions (const Json::Value &json_orbitals, const boost::shared_ptr<const Lattice> &lattice)
{
    const char * const json_orbitals_required[] = { "definitions", NULL };
    ensure_required(json_orbitals, json_orbitals_required);
    ensure_only(json_orbitals, json_orbitals_required);

    const Json::Value &json_defs = json_orbitals["definitions"];
    ensure_array(json_defs);

    const unsigned int N_filled = json_defs.size();
    const unsigned int N_sites = lattice->total_sites();

    if (N_filled > N_sites)
        throw ParseError("cannot have more orbitals than the number of sites on the lattice");

    Eigen::Matrix<amplitude_t, Eigen::Dynamic, Eigen::Dynamic> orbitals(N_filled, N_sites);

    for (unsigned int i = 0; i < N_filled; ++i) {
        const Json::Value &json_current_def = json_defs[i];
        ensure_array(json_current_def, N_sites);
        for (unsigned int j = 0; j < N_sites; ++j) {
            orbitals(i, j) = parse_complex(json_current_def[j]);
        }
    }

    return boost::make_shared<OrbitalDefinitions>(orbitals, lattice);
}

static boost::shared_ptr<const OrbitalDefinitions> parse_json_orbitals (const Json::Value &json_orbitals, const boost::shared_ptr<const Lattice > &lattice)
{
    return parse_json_orbitals_from_definitions(json_orbitals, lattice);
}

static boost::shared_ptr<const Subsystem> parse_json_subsystem (const Json::Value &json_subsystem, const Lattice &lattice)
{
    ensure_object_with_type_field_as_string(json_subsystem);
    if (std::strcmp(json_subsystem["type"].asCString(), "simple") == 0) {
        const char * const json_subsystem_required[] = { "type", "dimensions", NULL };
        ensure_required(json_subsystem, json_subsystem_required);
        ensure_only(json_subsystem, json_subsystem_required);
        const Json::Value &json_lengths = json_subsystem["dimensions"];
        const unsigned int n_dimensions = lattice.n_dimensions();
        ensure_array(json_lengths, n_dimensions);
        lw_vector<unsigned int, MAX_DIMENSION> lengths(n_dimensions);
        for (unsigned int i = 0; i < n_dimensions; ++i) {
            if (!(json_lengths[i].isIntegral() && json_lengths[i].asInt() >= 0))
                throw ParseError("subsystem length must be a non-negative integer");
            if (json_lengths[i].asInt() > lattice.dimensions[i])
                throw ParseError("subsystem length must fit within the lattice");
            lengths[i] = json_lengths[i].asUInt();
        }
        return boost::make_shared<SimpleSubsystem>(lengths);
    } else {
        throw ParseError("invalid subsystem type");
    }
}

static std::auto_ptr<RandomNumberGenerator> create_rng (const Json::Value &json_rng)
{
    unsigned long seed;
    ensure_object(json_rng);
    const char * const json_rng_allowed[] = { "seed", "type", NULL };
    ensure_only(json_rng, json_rng_allowed);
    if (json_rng.isMember("seed")) {
        if (!json_rng["seed"].isIntegral())
            throw ParseError("seed must be correct data type");
        seed = json_rng["seed"].asUInt();
    } else {
        throw ParseError("seed must be given");
    }
    const char *rng_type_name = "boost::mt19937"; // by default
    if (json_rng.isMember("type")) {
        const Json::Value &json_rng_type_name = json_rng["type"];
        ensure_string(json_rng_type_name);
        rng_type_name = json_rng_type_name.asCString();
        if (!RandomNumberGenerator::name_is_valid(rng_type_name))
            throw ParseError("invalid random number generator type specified");
    }
    return RandomNumberGenerator::create(rng_type_name, seed);
}

static boost::shared_ptr<Wavefunction> create_wavefunction (const Json::Value &json_wavefunction, const boost::shared_ptr<const Lattice> &lattice)
{
    ensure_object(json_wavefunction);
    const char * const json_wavefunction_required[] = { "type", NULL };
    ensure_required(json_wavefunction, json_wavefunction_required);
    ensure_string(json_wavefunction["type"]);
    const char *json_wavefunction_type_cstr = json_wavefunction["type"].asCString();
    if (std::strcmp(json_wavefunction_type_cstr, "free-fermion") == 0) {
        // free fermion wavefunction
        const char * const json_free_fermion_wavefunction_required[] = { "type", "orbitals", NULL };
        ensure_required(json_wavefunction, json_free_fermion_wavefunction_required);
        ensure_only(json_wavefunction, json_free_fermion_wavefunction_required);
        boost::shared_ptr<const OrbitalDefinitions> orbitals = parse_json_orbitals(json_wavefunction["orbitals"], lattice);
        return boost::make_shared<FreeFermionWavefunction>(orbitals);
    } else if (std::strcmp(json_wavefunction_type_cstr, "dbl") == 0) {
        // dbl wavefunction
        const char * const json_dbl_wavefunction_required[] = { "type", "orbitals-d1", "orbitals-d2", NULL };
        const char * const json_dbl_wavefunction_allowed[] = { "type", "orbitals-d1", "orbitals-d2", "exponent-d1", "exponent-d2", NULL };
        ensure_required(json_wavefunction, json_dbl_wavefunction_required);
        ensure_only(json_wavefunction, json_dbl_wavefunction_allowed);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d1 = parse_json_orbitals(json_wavefunction["orbitals-d1"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d2 = parse_json_orbitals(json_wavefunction["orbitals-d2"], lattice);
        if (orbitals_d1->get_N_filled() != orbitals_d2->get_N_filled())
            throw ParseError("d1 and d2 have different number of orbitals");
        return boost::make_shared<DBLWavefunction>(orbitals_d1, orbitals_d2,
                                                   json_get_double(json_wavefunction, "exponent-d1", 1.0),
                                                   json_get_double(json_wavefunction, "exponent-d2", 1.0));
    } else if (std::strcmp(json_wavefunction_type_cstr, "dmetal") == 0) {
        // dmetal wavefunction
        const char * const json_dmetal_wavefunction_required[] = { "type", "orbitals-d1", "orbitals-d2", "orbitals-f_up", "orbitals-f_down", NULL };
        const char * const json_dmetal_wavefunction_allowed[] = { "type", "orbitals-d1", "orbitals-d2", "orbitals-f_up", "orbitals-f_down", "exponent-d1", "exponent-d2", "exponent-f_up", "exponent-f_down", NULL };
        ensure_required(json_wavefunction, json_dmetal_wavefunction_required);
        ensure_only(json_wavefunction, json_dmetal_wavefunction_allowed);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d1 = parse_json_orbitals(json_wavefunction["orbitals-d1"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d2 = parse_json_orbitals(json_wavefunction["orbitals-d2"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_f_up = parse_json_orbitals(json_wavefunction["orbitals-f_up"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_f_down = parse_json_orbitals(json_wavefunction["orbitals-f_down"], lattice);
        if (orbitals_d1->get_N_filled() != orbitals_d2->get_N_filled())
            throw ParseError("d1 and d2 have different number of orbitals");
        if (orbitals_f_up->get_N_filled() + orbitals_f_down->get_N_filled() != orbitals_d1->get_N_filled())
            throw ParseError("number of orbitals in f_up and f_down must sum to number of orbitals in d1");
        return boost::make_shared<DMetalWavefunction>(orbitals_d1, orbitals_d2, orbitals_f_up, orbitals_f_down,
                                                      json_get_double(json_wavefunction, "exponent-d1", 1.0),
                                                      json_get_double(json_wavefunction, "exponent-d2", 1.0),
                                                      json_get_double(json_wavefunction, "exponent-f_up", 1.0),
                                                      json_get_double(json_wavefunction, "exponent-f_down", 1.0));
    } else if (std::strcmp(json_wavefunction_type_cstr, "rvb") == 0) {
        // rvb wavefunction
        const char * const json_rvb_wavefunction_required[] = { "type", "phi", NULL };
        ensure_required(json_wavefunction, json_rvb_wavefunction_required);
        ensure_only(json_wavefunction, json_rvb_wavefunction_required);
        const Json::Value &json_phi = json_wavefunction["phi"];
        const unsigned int N_sites = lattice->total_sites();
        if (N_sites % 2 == 1)
            throw ParseError("RVB wavefunction must be on a lattice with an even number of sites");
        ensure_array(json_phi, N_sites);
        std::vector<complex_t> phi(N_sites);
        for (unsigned int i = 0; i < N_sites; ++i)
            phi[i] = parse_complex(json_phi[i]);
        return boost::make_shared<RVBWavefunction>(lattice, phi);
    } else {
        throw ParseError("invalid wavefunction type");
    }
}

static std::auto_ptr<Walk> create_walk (const Json::Value &json_simulation, boost::shared_ptr<Wavefunction::Amplitude> &wfa)
{
    // begin setting up the simulation
    ensure_object(json_simulation);
    const char * const json_simulation_required[] = { "walk-type", NULL };
    ensure_required(json_simulation, json_simulation_required);

    // from here forward, we have special logic per walk type
    ensure_string(json_simulation["walk-type"]);
    const char *json_walk_type_cstr = json_simulation["walk-type"].asCString();
    if (std::strcmp(json_walk_type_cstr, "standard") == 0) {

        // STANDARD WALK

        // ensure correct json properties are given
        const char * const json_simulation_only[] = { "walk-type", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up walk
        return std::auto_ptr<Walk>(new StandardWalk(wfa));

    } else if (std::strcmp(json_walk_type_cstr, "renyi-mod/possible") == 0) {

        // RENYI MOD/POSSIBLE WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { "walk-type", "subsystem", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem(json_simulation["subsystem"], wfa->get_lattice()));

        // We need two copies of the system, each of which has the same number
        // of particles in the subsystem.  So for now we just initialize both
        // copies with the same exact positions.
        boost::shared_ptr<Wavefunction::Amplitude> wfa2(wfa->clone());

        // set up walk
        return std::auto_ptr<Walk>(new RenyiModPossibleWalk(wfa, wfa2, subsystem));

    } else if (std::strcmp(json_walk_type_cstr, "renyi-sign") == 0) {

        // RENYI SIGN WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { "walk-type", "subsystem", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem(json_simulation["subsystem"], wfa->get_lattice()));

        // We need two copies of the system, each of which has the same number
        // of particles in the subsystem.  So for now we just initialize both
        // copies with the same exact positions.
        boost::shared_ptr<Wavefunction::Amplitude> wfa2(wfa->clone());

        // set up walk
        return std::auto_ptr<Walk>(new RenyiSignWalk(wfa, wfa2, subsystem));

    } else {
        throw ParseError("invalid walk type");
    }
}

std::auto_ptr<MetropolisSimulation> create_simulation (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice, const std::list<boost::shared_ptr<BaseMeasurement> > &measurements, unsigned int equilibrium_steps)
{
    Json::Value json_input;
    {
        std::istringstream iss(json_input_str, std::istringstream::in);
        iss >> json_input;
    }

    ensure_object(json_input);
    const char * const json_input_required[] = { "rng", "system", "simulation", NULL };
    ensure_required(json_input, json_input_required);
    ensure_only(json_input, json_input_required);

    // initialize random number generator
    std::auto_ptr<RandomNumberGenerator> rng(create_rng(json_input["rng"]));

    // set up the wavefunction
    const Json::Value &json_system = json_input["system"];
    ensure_object(json_system);
    const char * const json_system_required[] = { "wavefunction", NULL };
    ensure_required(json_system, json_system_required);
    ensure_only(json_system, json_system_required);
    boost::shared_ptr<Wavefunction> wf(create_wavefunction(json_system["wavefunction"], lattice));

    // set up the Wavefunction::Amplitude
    boost::shared_ptr<Wavefunction::Amplitude> wfa(wf->create_nonzero_wavefunctionamplitude(wf, *rng));
    if (!wfa)
        throw ParseError("could not find a nonzero wavefunction amplitude");

    // set up the walk
    std::auto_ptr<Walk> walk(create_walk(json_input["simulation"], wfa));

    // return the simulation
    return std::auto_ptr<MetropolisSimulation>(new MetropolisSimulation(walk, measurements, equilibrium_steps, rng));
}
