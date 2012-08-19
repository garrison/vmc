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
#include "random-configuration.hpp"
#include "PositionArguments.hpp"
#include "OrbitalDefinitions.hpp"
#include "FreeFermionWavefunctionAmplitude.hpp"
#include "DBLWavefunctionAmplitude.hpp"
#include "DMetalWavefunctionAmplitude.hpp"
#include "RVBWavefunctionAmplitude.hpp"
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
#include "RunInformation.hpp"
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

static unsigned int parse_uint (const Json::Value &uint)
{
    if (!uint.isIntegral() || uint.asInt() < 0)
        throw ParseError("unsigned integer expected");
    return uint.asUInt();
}

static unsigned int parse_uint (const Json::Value &uint, unsigned int lowest_out_of_range)
{
    const unsigned int rv = parse_uint(uint);
    // fixme: this error message is not very descriptive ...
    if (rv >= lowest_out_of_range)
        throw ParseError("expecting an unsigned integer in a certain range");
    return rv;
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

static unsigned int parse_json_steps_per_measurement (const Json::Value &json_measurement_def)
{
    if (json_measurement_def.isMember("steps-per-measurement")) {
        const Json::Value &spm = json_measurement_def["steps-per-measurement"];
        if (!(spm.isIntegral() && spm.asInt() > 0))
            throw ParseError("invalid steps-per-measurement value");
        return spm.asUInt();
    } else {
        return 1;
    }
}

static void ensure_single_step_per_measurement (const Json::Value &json_measurement_def)
{
    if (parse_json_steps_per_measurement(json_measurement_def) != 1)
        throw ParseError("given measurement requires a measurement on each step");
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

static BoundaryConditions parse_json_boundary_conditions (const Json::Value &json_bcs, unsigned int n_dimensions)
{
    ensure_array(json_bcs, n_dimensions);
    BoundaryConditions boundary_conditions(n_dimensions);
    for (unsigned int i = 0; i < n_dimensions; ++i) {
        if (!(json_bcs[i].isIntegral() && json_bcs[i].asInt() > 0))
            throw ParseError("invalid boundary condition specifier");
        boundary_conditions[i] = BoundaryCondition(json_bcs[i].asUInt());
    }
    return boundary_conditions;
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

static boost::shared_ptr<Measurement<StandardWalk> > parse_standard_walk_measurement_definition (const Json::Value &json_measurement_def, const WavefunctionAmplitude &wf, const boost::shared_ptr<const Lattice> &lattice)
{
    const unsigned int N_species = wf.get_positions().get_N_species();
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "operator") == 0) {
        const char * const json_operator_required[] = { "type", "hops", NULL };
        const char * const json_operator_allowed[] = { "type", "hops", "sum", "boundary-conditions", "steps-per-measurement", NULL };
        ensure_required(json_measurement_def, json_operator_required);
        ensure_only(json_measurement_def, json_operator_allowed);
        unsigned int steps_per_measurement = parse_json_steps_per_measurement(json_measurement_def);
        const Json::Value &json_hops = json_measurement_def["hops"];
        ensure_array(json_hops);
        std::vector<SiteHop> site_hop_v;
        for (unsigned int i = 0; i < json_hops.size(); ++i) {
            const Json::Value &json_current_hop = json_hops[i];
            ensure_object(json_current_hop);
            const char * const json_current_hop_required[] = { "source", "destination", "species", NULL };
            ensure_required(json_current_hop, json_current_hop_required);
            ensure_only(json_current_hop, json_current_hop_required);
            const unsigned int N_sites = wf.get_positions().get_N_sites();
            const unsigned int N_species = wf.get_positions().get_N_species();
            const LatticeSite source(lattice->site_from_index(parse_uint(json_current_hop["source"], N_sites)));
            const LatticeSite destination(lattice->site_from_index(parse_uint(json_current_hop["destination"], N_sites)));
            const unsigned int species = parse_uint(json_current_hop["species"], N_species);
            site_hop_v.push_back(SiteHop(source, destination, species));
        }
        if (!ParticleOperator::is_valid(site_hop_v, *lattice, N_species))
            throw ParseError("invalid operator provided for operator measurement");
        bool sum = false;
        if (json_measurement_def.isMember("sum")) {
            if (!json_measurement_def["sum"].isBool())
                throw ParseError("boolean expected");
            sum = json_measurement_def["sum"].asBool();
        }
        boost::shared_ptr<const BoundaryConditions> bcs;
        if (json_measurement_def.isMember("boundary-conditions") && !json_measurement_def["boundary-conditions"].isNull()) {
            if (!sum)
                throw ParseError("if boundary conditions are given for an operator, it must be a sum over sites");
            bcs.reset(new BoundaryConditions(parse_json_boundary_conditions(json_measurement_def["boundary-conditions"],
                                                                            lattice->n_dimensions())));
        }
        return boost::make_shared<OperatorMeasurement>(steps_per_measurement, ParticleOperator(site_hop_v, lattice), sum, bcs.get());
    } else if (std::strcmp(json_measurement_def["type"].asCString(), "subsystem-occupation-number-probability") == 0) {
        const char * const json_subsystem_occupation_required[] = { "type", "subsystem", NULL };
        const char * const json_subsystem_occupation_allowed[] = { "type", "subsystem", "steps-per-measurement", NULL };
        ensure_required(json_measurement_def, json_subsystem_occupation_required);
        ensure_only(json_measurement_def, json_subsystem_occupation_allowed);
        unsigned int steps_per_measurement = parse_json_steps_per_measurement(json_measurement_def);
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem(json_measurement_def["subsystem"], *lattice));
        return boost::make_shared<SubsystemOccupationNumberProbabilityMeasurement>(steps_per_measurement, subsystem);
    } else {
        throw ParseError("invalid standard walk measurement type");
    }
}

static boost::shared_ptr<Measurement<RenyiModPossibleWalk> > parse_renyi_mod_possible_walk_measurement_definition (const Json::Value &json_measurement_def)
{
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "renyi-mod/possible") == 0) {
        const char * const json_renyi_mod_possible_allowed[] = { "type", "steps-per-measurement", NULL };
        ensure_only(json_measurement_def, json_renyi_mod_possible_allowed);
        ensure_single_step_per_measurement(json_measurement_def);
        return boost::make_shared<RenyiModPossibleMeasurement>();
    } else {
        throw ParseError("invalid renyi-mod/possible walk measurement type");
    }
}

static boost::shared_ptr<Measurement<RenyiSignWalk> > parse_renyi_sign_walk_measurement_definition (const Json::Value &json_measurement_def)
{
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "renyi-sign") == 0) {
        const char * const json_renyi_sign_allowed[] = { "type", "steps-per-measurement", NULL };
        ensure_only(json_measurement_def, json_renyi_sign_allowed);
        ensure_single_step_per_measurement(json_measurement_def);
        return boost::make_shared<RenyiSignMeasurement>();
    } else {
        throw ParseError("invalid renyi-sign walk measurement type");
    }
}

static void set_wavefunction_positions_from_json (WavefunctionAmplitude &wf, const Json::Value &json_initial_positions)
{
    const unsigned int N_sites = wf.get_positions().get_N_sites();
    const unsigned int N_species = wf.get_positions().get_N_species();
    ensure_array(json_initial_positions, N_species);
    std::vector<std::vector<unsigned int> > v;
    for (unsigned int species = 0; species < N_species; ++species) {
        const unsigned int N_filled = wf.get_positions().get_N_filled(species);
        ensure_array(json_initial_positions[species], N_filled);
        std::set<unsigned int> vs;
        v.push_back(std::vector<unsigned int>(N_filled));
        for (unsigned int i = 0; i < N_filled; ++i) {
            const Json::Value &json_pos = json_initial_positions[species][i];
            if (!json_pos.isIntegral())
                throw ParseError("expecting integer");
            if (!(json_pos.asInt() >= 0 && json_pos.asUInt() < N_sites))
                throw ParseError("invalid site index");
            bool inserted = vs.insert(json_pos.asUInt()).second;
            if (!inserted)
                throw ParseError("position index specified twice, but double occupancy is not allowed");
            v[species][i] = json_pos.asUInt();
        }
    }
    wf.reset(PositionArguments(v, N_sites));
    if (wf.psi() == amplitude_t(0))
        throw ParseError("given positions have zero amplitude");
}

static Json::Value positions_json_repr (const PositionArguments &r)
{
    Json::Value rv(Json::arrayValue);
    for (unsigned int species = 0; species < r.get_N_species(); ++species) {
        Json::Value sv(Json::arrayValue);
        for (unsigned int i = 0; i < r.get_N_filled(species); ++i) {
            sv.append(r[Particle(i, species)]);
        }
        rv.append(sv);
    }
    return rv;
}

static Json::Value monte_carlo_stats_json_repr (const MetropolisSimulation &sim)
{
    Json::Value rv(Json::objectValue);
    rv["steps-completed"] = Json::UInt(sim.steps_completed());
    rv["steps-accepted"] = Json::UInt(sim.steps_accepted());
    rv["steps-fully-rejected"] = Json::UInt(sim.steps_fully_rejected());
    return rv;
}

static inline double jsoncpp_real_cast (real_t v)
{
    // it would be really nice if jsoncpp supported "long double" directly...
    // (we use this cast explicitly so that if we ever have real support for
    // long double, we can add it by grepping for calls of this function)
    return v;
}

static Json::Value complex_to_json_array (const complex_t &v)
{
    Json::Value rv(Json::arrayValue);
    rv.append(Json::Value(jsoncpp_real_cast(std::real(v))));
    rv.append(Json::Value(jsoncpp_real_cast(std::imag(v))));
    return rv;
}

static Json::Value standard_walk_measurement_json_repr (const BaseMeasurement *measurement_ptr, const WavefunctionAmplitude &wf)
{
    const OperatorMeasurement *om = dynamic_cast<const OperatorMeasurement *>(measurement_ptr);
    const SubsystemOccupationNumberProbabilityMeasurement *sonpm = dynamic_cast<const SubsystemOccupationNumberProbabilityMeasurement *>(measurement_ptr);
    if (om) {
        // operator measurement
        return Json::Value(complex_to_json_array(om->get()));
    } else if (sonpm) {
        // subsystem occupation number probability measurement
        const PositionArguments &r = wf.get_positions();
        Json::Value rv(Json::arrayValue);
        std::vector<unsigned int> occupation(r.get_N_species());
        Json::Value current_pair(Json::arrayValue);
        current_pair.resize(2);
        current_pair[0].resize(occupation.size());
        bool done = false;
        while (!done) {
            // append to json the current occupation probability (if nonzero)
            real_t val = sonpm->get(occupation);
            if (val != real_t(0)) {
                for (unsigned int i = 0; i < occupation.size(); ++i)
                    current_pair[0][i] = occupation[i];
                current_pair[1] = jsoncpp_real_cast(val);
                rv.append(current_pair);
            }

            // advance to the next occupation we should consider
            for (unsigned int j = 0; j < occupation.size(); ++j) {
                ++occupation[j];
                if (occupation[j] == r.get_N_filled(j)) {
                    occupation[j] = 0;
                    if (j == occupation.size() - 1)
                        done = true;
                } else {
                    break;
                }
            }
        }
        return rv;
    } else {
        // should not be reached
        BOOST_ASSERT(false);
        return Json::Value();
    }
}

static Json::Value renyi_mod_possible_walk_measurement_json_repr (const BaseMeasurement *measurement_ptr)
{
    const RenyiModPossibleMeasurement *rmpm = boost::polymorphic_downcast<const RenyiModPossibleMeasurement*>(measurement_ptr);
    return Json::Value(jsoncpp_real_cast(rmpm->get()));
}

static Json::Value renyi_sign_walk_measurement_json_repr (const BaseMeasurement *measurement_ptr)
{
    const RenyiSignMeasurement *rsm = boost::polymorphic_downcast<const RenyiSignMeasurement*>(measurement_ptr);
    return complex_to_json_array(rsm->get());
}

HighlevelSimulation::HighlevelSimulation (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice)
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
    unsigned long seed;
    const Json::Value &json_rng = json_input["rng"];
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
    std::auto_ptr<RandomNumberGenerator> rng(RandomNumberGenerator::create(rng_type_name, seed));

    // begin setting up the physical system
    const Json::Value &json_system = json_input["system"];
    ensure_object(json_system);
    const char * const json_system_required[] = { "wavefunction", NULL };
    ensure_required(json_system, json_system_required);
    ensure_only(json_system, json_system_required);

    // set up the wavefunction
    boost::shared_ptr<WavefunctionAmplitude> wf;
    const Json::Value &json_wavefunction = json_input["system"]["wavefunction"];
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
        std::vector<std::vector<unsigned int> > configuration;
        configuration.push_back(some_random_configuration(orbitals->get_N_filled(), *lattice, *rng));
        wf.reset(new FreeFermionWavefunctionAmplitude(PositionArguments(configuration, lattice->total_sites()), orbitals));
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
        std::vector<std::vector<unsigned int> > configuration;
        configuration.push_back(some_random_configuration(orbitals_d1->get_N_filled(), *lattice, *rng));
        wf.reset(new DBLWavefunctionAmplitude(PositionArguments(configuration, lattice->total_sites()),
                                              orbitals_d1, orbitals_d2,
                                              json_get_double(json_wavefunction, "exponent-d1", 1.0),
                                              json_get_double(json_wavefunction, "exponent-d2", 1.0)));
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
        std::vector<std::vector<unsigned int> > configuration;
        configuration.push_back(some_random_configuration(orbitals_f_up->get_N_filled(), *lattice, *rng));
        configuration.push_back(some_random_configuration(orbitals_f_down->get_N_filled(), *lattice, *rng));
        wf.reset(new DMetalWavefunctionAmplitude(PositionArguments(configuration, lattice->total_sites()),
                                                 orbitals_d1, orbitals_d2, orbitals_f_up, orbitals_f_down,
                                                 json_get_double(json_wavefunction, "exponent-d1", 1.0),
                                                 json_get_double(json_wavefunction, "exponent-d2", 1.0),
                                                 json_get_double(json_wavefunction, "exponent-f_up", 1.0),
                                                 json_get_double(json_wavefunction, "exponent-f_down", 1.0)));
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

        // find some initial positions with no double occupancy
        std::vector<std::vector<unsigned int> > configuration;
        configuration.push_back(some_random_configuration(N_sites / 2, *lattice, *rng));
        configuration.push_back(std::vector<unsigned int>());
        std::vector<unsigned int> occupied_sites(N_sites);
        for (unsigned int i = 0; i < configuration[0].size(); ++i) {
            BOOST_ASSERT(configuration[0][i] < N_sites);
            BOOST_ASSERT(occupied_sites[configuration[0][i]] == 0);
            ++occupied_sites[configuration[0][i]];
        }
        for (unsigned int i = 0; i < N_sites; ++i) {
            if (occupied_sites[i] == 0)
                configuration[1].push_back(i);
        }
        BOOST_ASSERT(configuration[0].size() == configuration[1].size());
        wf.reset(new RVBWavefunctionAmplitude(PositionArguments(configuration, lattice->total_sites()), lattice, phi));
    } else {
        throw ParseError("invalid wavefunction type");
    }

    // begin setting up the simulation
    const Json::Value &json_simulation = json_input["simulation"];
    ensure_object(json_simulation);
    const char * const json_simulation_required[] = { "walk-type", NULL };
    ensure_required(json_simulation, json_simulation_required);
#define JSON_SIMULATION_GLOBAL_ALLOWED "walk-type", "equilibrium-steps", "initial-positions"

    // determine how many steps to do
    unsigned int equilibrium_steps = 0;
    if (json_simulation.isMember("equilibrium-steps")) {
        const Json::Value &json_equilibrium_steps = json_simulation["equilibrium-steps"];
        if (!(json_equilibrium_steps.isIntegral() && json_equilibrium_steps.asInt() >= 0))
            throw ParseError("equilibrium-steps must be a non-negative integer");
        equilibrium_steps = json_equilibrium_steps.asUInt();
    }

    // if there are no equilibrium steps, make sure initial positions are given
    if (equilibrium_steps == 0 && !json_simulation.isMember("initial-positions"))
        throw ParseError("initial positions must be given if there are no equilibrium steps");

    // from here forward, we have special logic per walk type
    ensure_string(json_simulation["walk-type"]);
    const char *json_walk_type_cstr = json_simulation["walk-type"].asCString();
    if (std::strcmp(json_walk_type_cstr, "standard") == 0) {

        // STANDARD WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "measurements", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "measurements", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up initial positions of wavefunction
        if (json_simulation.isMember("initial-positions")) {
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"]);
            if (wf->psi() == amplitude_t(0))
                throw ParseError("wavefunction has zero amplitude at given initial-positions");
        } else {
            bool success = search_for_configuration_with_nonzero_amplitude(*wf, *rng);
            if (!success)
                throw ParseError("could not find a configuration with nonzero amplitude");
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<BaseMeasurement> > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i) {
            measurements.push_back(parse_standard_walk_measurement_definition(json_measurements[i], *wf, lattice));
        }

        // set up and perform walk
        std::auto_ptr<Walk> walk(new StandardWalk(wf));
        sim.reset(new MetropolisSimulation(walk, measurements, equilibrium_steps, rng));

    } else if (std::strcmp(json_walk_type_cstr, "renyi-mod/possible") == 0) {

        // RENYI MOD/POSSIBLE WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", "measurements", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "subsystem", "measurements", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem(json_simulation["subsystem"], *lattice));

        // set up initial positions of wavefunctions
        boost::shared_ptr<WavefunctionAmplitude> wf2(wf->clone());
        if (json_simulation.isMember("initial-positions")) {
            ensure_array(json_simulation["initial-positions"], 2);
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"][0u]);
            set_wavefunction_positions_from_json(*wf2, json_simulation["initial-positions"][1u]);
            if (!count_subsystem_particle_counts_for_match(*wf, *wf2, *subsystem))
                throw ParseError("The initial positions of each copy must have the same numbers/types of particles in the subsystem");
        } else {
            bool success = search_for_configuration_with_nonzero_amplitude(*wf, *rng);
            if (!success)
                throw ParseError("could not find a configuration with nonzero amplitude");
            // We need two copies of the system, each of which has the same
            // number of particles in the subsystem.  So for now we just
            // initialize both copies with the same exact positions.
            wf2 = wf;
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<BaseMeasurement> > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i) {
            measurements.push_back(parse_renyi_mod_possible_walk_measurement_definition(json_measurements[i]));
        }

        // set up and perform walk
        std::auto_ptr<Walk> walk(new RenyiModPossibleWalk(wf, wf2, subsystem));
        sim.reset(new MetropolisSimulation(walk, measurements, equilibrium_steps, rng));

    } else if (std::strcmp(json_walk_type_cstr, "renyi-sign") == 0) {

        // RENYI SIGN WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", "measurements", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "subsystem", "measurements", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem(json_simulation["subsystem"], *lattice));

        // set up initial positions of wavefunctions
        boost::shared_ptr<WavefunctionAmplitude> wf2(wf->clone());
        if (json_simulation.isMember("initial-positions")) {
            ensure_array(json_simulation["initial-positions"], 2);
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"][0u]);
            set_wavefunction_positions_from_json(*wf2, json_simulation["initial-positions"][1u]);
            if (!count_subsystem_particle_counts_for_match(*wf, *wf2, *subsystem))
                throw ParseError("The initial positions of each copy must have the same numbers/types of particles in the subsystem");
        } else {
            bool success = search_for_configuration_with_nonzero_amplitude(*wf, *rng);
            if (!success)
                throw ParseError("could not find a configuration with nonzero amplitude");
            // We need two copies of the system, each of which has the same
            // number of particles in the subsystem.  So for now we just
            // initialize both copies with the same exact positions.
            wf2 = wf;
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<BaseMeasurement> > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i) {
            measurements.push_back(parse_renyi_sign_walk_measurement_definition(json_measurements[i]));
        }

        // set up and perform walk
        std::auto_ptr<Walk> walk(new RenyiSignWalk(wf, wf2, subsystem));
        sim.reset(new MetropolisSimulation(walk, measurements, equilibrium_steps, rng));

    } else {
        throw ParseError("invalid walk type");
    }

    walk_type = json_walk_type_cstr;
}

std::string HighlevelSimulation::output (void) const
{
    Json::Value json_measurement_output(Json::arrayValue);
    Json::Value json_final_positions_output, json_monte_carlo_stats_output;
    const std::list<boost::shared_ptr<BaseMeasurement> > &measurements = sim->get_measurements();
    if (walk_type == "standard") {
        const StandardWalk *walk_ptr = static_cast<const StandardWalk *>(sim->get_walk_ptr());
        for (std::list<boost::shared_ptr<BaseMeasurement> >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(standard_walk_measurement_json_repr(i->get(), walk_ptr->get_wavefunction()));
        }
        json_final_positions_output = positions_json_repr(walk_ptr->get_wavefunction().get_positions());
    } else if (walk_type == "renyi-mod/possible") {
        const RenyiModPossibleWalk *walk_ptr = static_cast<const RenyiModPossibleWalk *>(sim->get_walk_ptr());
        for (std::list<boost::shared_ptr<BaseMeasurement> >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(renyi_mod_possible_walk_measurement_json_repr(i->get()));
        }
        json_final_positions_output = Json::Value(Json::arrayValue);
        json_final_positions_output.append(positions_json_repr(walk_ptr->get_phialpha1().get_positions()));
        json_final_positions_output.append(positions_json_repr(walk_ptr->get_phialpha2().get_positions()));
    } else if (walk_type == "renyi-sign") {
        const RenyiSignWalk *walk_ptr = static_cast<const RenyiSignWalk *>(sim->get_walk_ptr());
        for (std::list<boost::shared_ptr<BaseMeasurement> >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(renyi_sign_walk_measurement_json_repr(i->get()));
        }
        json_final_positions_output = Json::Value(Json::arrayValue);
        json_final_positions_output.append(positions_json_repr(walk_ptr->get_phialpha1().get_positions()));
        json_final_positions_output.append(positions_json_repr(walk_ptr->get_phialpha2().get_positions()));
    } else {
        BOOST_ASSERT(false);
    }
    json_monte_carlo_stats_output = monte_carlo_stats_json_repr(*sim);

    Json::Value json_output(Json::objectValue);
    json_output["final-positions"] = json_final_positions_output;
    json_output["measurements"] = json_measurement_output;
    json_output["run-information"] = RunInformation::json_info();
    json_output["monte-carlo-stats"] = json_monte_carlo_stats_output;

    std::ostringstream oss(std::ostringstream::out);
    oss << json_output;
    return oss.str();
}
