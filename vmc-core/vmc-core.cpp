/*
 * vmc-core.cpp: everything necessary to set up a simulation based on json
 * input, run the simulation, and return its results in json format.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
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
#include <Eigen/Dense>

#include "NDLattice.hpp"
#include "random-filling.hpp"
#include "PositionArguments.hpp"
#include "OrbitalDefinitions.hpp"
#include "FilledOrbitals.hpp"
#include "FreeFermionWavefunctionAmplitude.hpp"
#include "DBLWavefunctionAmplitude.hpp"
#include "DMetalWavefunctionAmplitude.hpp"
#include "RVBWavefunctionAmplitude.hpp"
#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "DensityDensityMeasurement.hpp"
#include "GreenMeasurement.hpp"
#include "SubsystemOccupationNumberProbabilityMeasurement.hpp"
#include "RenyiModWalk.hpp"
#include "RenyiModMeasurement.hpp"
#include "RenyiModPossibleWalk.hpp"
#include "RenyiModPossibleMeasurement.hpp"
#include "RenyiSignWalk.hpp"
#include "RenyiSignMeasurement.hpp"
#include "Subsystem.hpp"
#include "SimpleSubsystem.hpp"
#include "RunInformation.hpp"

static RunInformation run_information;

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

template <unsigned int DIM>
boost::shared_ptr<const OrbitalDefinitions> parse_json_orbitals_from_definitions (const Json::Value &json_orbitals, const boost::shared_ptr<const NDLattice<DIM> > &lattice)
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

template <unsigned int DIM>
boost::shared_ptr<const OrbitalDefinitions> parse_json_orbitals_from_filling (const Json::Value &json_orbitals, const boost::shared_ptr<const NDLattice<DIM> > &lattice)
{
    const char * const json_orbitals_required[] = { "filling", "boundary-conditions", NULL };
    ensure_required(json_orbitals, json_orbitals_required);
    ensure_only(json_orbitals, json_orbitals_required);

    // set up the boundary conditions
    const Json::Value &json_bcs = json_orbitals["boundary-conditions"];
    ensure_array(json_bcs, DIM);
    typename NDLattice<DIM>::BoundaryConditions boundary_conditions;
    for (unsigned int i = 0; i < DIM; ++i) {
        if (!(json_bcs[i].isIntegral() && json_bcs[i].asInt() > 0))
            throw ParseError("invalid boundary condition specifier");
        boundary_conditions[i] = BoundaryCondition(json_bcs[i].asUInt());
    }

    // set up the orbitals' filled momenta
    const Json::Value &json_filling = json_orbitals["filling"];
    ensure_array(json_filling);
    std::vector<boost::array<int, DIM> > filled_momenta;
    std::set<boost::array<int, DIM> > filled_momenta_set;
    filled_momenta.reserve(json_filling.size());
    for (unsigned int i = 0; i < json_filling.size(); ++i) {
        const Json::Value &json_current_filling = json_filling[i];
        ensure_array(json_current_filling, DIM);
        boost::array<int, DIM> current_filling;
        for (unsigned int j = 0; j < DIM; ++j) {
            if (!(json_current_filling[j].isIntegral() && json_current_filling[j].asInt() >= 0 && json_current_filling[j].asInt() < lattice->length[j]))
                throw ParseError("invalid momentum index");
            current_filling[j] = json_current_filling[j].asInt();
        }
        filled_momenta.push_back(current_filling);
        if (!filled_momenta_set.insert(current_filling).second)
            throw ParseError("momentum was specified twice for the same orbital definitions");
    }

    return boost::make_shared<FilledOrbitals<DIM> >(filled_momenta, lattice, boundary_conditions);
}

template <unsigned int DIM>
boost::shared_ptr<const OrbitalDefinitions> parse_json_orbitals (const Json::Value &json_orbitals, const boost::shared_ptr<const NDLattice<DIM> > &lattice)
{
    if (json_orbitals.isMember("definitions")) {
        return parse_json_orbitals_from_definitions<DIM>(json_orbitals, lattice);
    } else {
        return parse_json_orbitals_from_filling<DIM>(json_orbitals, lattice);
    }
}

template <unsigned int DIM>
boost::shared_ptr<const Subsystem> parse_json_subsystem (const Json::Value &json_subsystem, const NDLattice<DIM> &lattice)
{
    ensure_object_with_type_field_as_string(json_subsystem);
    if (std::strcmp(json_subsystem["type"].asCString(), "simple") == 0) {
        const char * const json_subsystem_required[] = { "type", "dimensions", NULL };
        ensure_required(json_subsystem, json_subsystem_required);
        ensure_only(json_subsystem, json_subsystem_required);
        const Json::Value &json_lengths = json_subsystem["dimensions"];
        ensure_array(json_lengths, DIM);
        boost::array<unsigned int, DIM> lengths;
        for (unsigned int i = 0; i < DIM; ++i) {
            if (!(json_lengths[i].isIntegral() && json_lengths[i].asInt() >= 0))
                throw ParseError("subsystem length must be a non-negative integer");
            if (json_lengths[i].asInt() > lattice.length[i])
                throw ParseError("subsystem length must fit within the lattice");
            lengths[i] = json_lengths[i].asUInt();
        }
        return boost::make_shared<SimpleSubsystem<DIM> >(lengths);
    } else {
        throw ParseError("invalid subsystem type");
    }
}

template <unsigned int DIM>
boost::shared_ptr<Measurement<StandardWalk> > parse_standard_walk_measurement_definition (const Json::Value &json_measurement_def, const WavefunctionAmplitude &wf, const NDLattice<DIM> &lattice)
{
    const unsigned int N_species = wf.get_positions().get_N_species();
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "density-density") == 0) {
        const char * const json_density_density_allowed[] = { "type", "steps-per-measurement", "species", NULL };
        ensure_only(json_measurement_def, json_density_density_allowed);
        unsigned int steps_per_measurement = parse_json_steps_per_measurement(json_measurement_def);
        unsigned int species1 = 0, species2 = 0;
        if (json_measurement_def.isMember("species")) {
            if (json_measurement_def["species"].isArray()) {
                ensure_array(json_measurement_def["species"], 2);
                species1 = parse_uint(json_measurement_def["species"][0], N_species);
                species2 = parse_uint(json_measurement_def["species"][1], N_species);
            } else {
                species1 = parse_uint(json_measurement_def["species"], N_species);
                species2 = species1;
            }
        } else if (N_species != 1) {
            throw ParseError("species must be given, as this wave function contains multiple species of particles");
        }
        return boost::make_shared<DensityDensityMeasurement<DIM> >(steps_per_measurement, species1, species2);
    } else if (std::strcmp(json_measurement_def["type"].asCString(), "green") == 0) {
        const char * const json_density_density_allowed[] = { "type", "steps-per-measurement", "species", NULL };
        ensure_only(json_measurement_def, json_density_density_allowed);
        unsigned int steps_per_measurement = parse_json_steps_per_measurement(json_measurement_def);
        unsigned int species = 0;
        if (json_measurement_def.isMember("species")) {
            species = parse_uint(json_measurement_def["species"], N_species);
        } else if (N_species != 1) {
            throw ParseError("species must be given, as this wave function contains multiple species of particles");
        }
        return boost::make_shared<GreenMeasurement<DIM> >(steps_per_measurement, species);
    } else if (std::strcmp(json_measurement_def["type"].asCString(), "subsystem-occupation-number-probability") == 0) {
        const char * const json_subsystem_occupation_required[] = { "type", "subsystem", NULL };
        const char * const json_subsystem_occupation_allowed[] = { "type", "subsystem", "steps-per-measurement", NULL };
        ensure_required(json_measurement_def, json_subsystem_occupation_required);
        ensure_only(json_measurement_def, json_subsystem_occupation_allowed);
        unsigned int steps_per_measurement = parse_json_steps_per_measurement(json_measurement_def);
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem<DIM>(json_measurement_def["subsystem"], lattice));
        return boost::make_shared<SubsystemOccupationNumberProbabilityMeasurement>(subsystem, steps_per_measurement);
    } else {
        throw ParseError("invalid standard walk measurement type");
    }
}

template <unsigned int DIM>
boost::shared_ptr<Measurement<RenyiModWalk> > parse_renyi_mod_walk_measurement_definition (const Json::Value &json_measurement_def, const NDLattice<DIM> &lattice)
{
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "renyi-mod") == 0) {
        const char * const json_renyi_mod_required[] = { "type", "subsystem", NULL };
        const char * const json_renyi_mod_allowed[] = { "type", "steps-per-measurement", "subsystem", NULL };
        ensure_required(json_measurement_def, json_renyi_mod_required);
        ensure_only(json_measurement_def, json_renyi_mod_allowed);
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem<DIM>(json_measurement_def["subsystem"], lattice));
        unsigned int steps_per_measurement = parse_json_steps_per_measurement(json_measurement_def);
        return boost::make_shared<RenyiModMeasurement>(subsystem, steps_per_measurement);
    } else {
        throw ParseError("invalid renyi-mod walk measurement type");
    }
}

template <unsigned int DIM>
boost::shared_ptr<Measurement<RenyiModPossibleWalk> > parse_renyi_mod_possible_walk_measurement_definition (const Json::Value &json_measurement_def)
{
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "renyi-mod/possible") == 0) {
        const char * const json_renyi_mod_possible_allowed[] = { "type", NULL };
        ensure_only(json_measurement_def, json_renyi_mod_possible_allowed);
        return boost::make_shared<RenyiModPossibleMeasurement>();
    } else {
        throw ParseError("invalid renyi-mod/possible walk measurement type");
    }
}

template <unsigned int DIM>
boost::shared_ptr<Measurement<RenyiSignWalk> > parse_renyi_sign_walk_measurement_definition (const Json::Value &json_measurement_def)
{
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "renyi-sign") == 0) {
        const char * const json_renyi_sign_allowed[] = { "type", NULL };
        ensure_only(json_measurement_def, json_renyi_sign_allowed);
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

template <class Walk_T>
static Json::Value monte_carlo_stats_json_repr (const MetropolisSimulation<Walk_T> &sim)
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

template <unsigned int DIM>
static Json::Value standard_walk_measurement_json_repr (const Measurement<StandardWalk> *measurement_ptr, const WavefunctionAmplitude &wf)
{
    const DensityDensityMeasurement<DIM> *ddm = dynamic_cast<const DensityDensityMeasurement<DIM>*>(measurement_ptr);
    const GreenMeasurement<DIM> *gm = dynamic_cast<const GreenMeasurement<DIM>*>(measurement_ptr);
    const SubsystemOccupationNumberProbabilityMeasurement *sonpm = dynamic_cast<const SubsystemOccupationNumberProbabilityMeasurement *>(measurement_ptr);
    if (ddm) {
        // density-density measurement
        Json::Value rv(Json::arrayValue);
        for (unsigned int i = 0; i < ddm->basis_indices(); ++i) {
            Json::Value a(Json::arrayValue);
            for (unsigned int j = 0; j < ddm->get_N_sites(); ++j)
                a.append(Json::Value(jsoncpp_real_cast(ddm->get(j, i))));
            rv.append(a);
        }
        return rv;
    } else if (gm) {
        // green measurement
        Json::Value rv(Json::arrayValue);
        for (unsigned int i = 0; i < gm->basis_indices(); ++i) {
            Json::Value a(Json::arrayValue);
            for (unsigned int j = 0; j < gm->get_N_sites(); ++j)
                a.append(Json::Value(complex_to_json_array(gm->get(j, i))));
            rv.append(a);
        }
        return rv;
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

static Json::Value renyi_mod_walk_measurement_json_repr (const Measurement<RenyiModWalk> *measurement_ptr)
{
    const RenyiModMeasurement *rmm = boost::polymorphic_downcast<const RenyiModMeasurement*>(measurement_ptr);
    return Json::Value(rmm->get());
}

static Json::Value renyi_mod_possible_walk_measurement_json_repr (const Measurement<RenyiModPossibleWalk> *measurement_ptr)
{
    const RenyiModPossibleMeasurement *rmpm = boost::polymorphic_downcast<const RenyiModPossibleMeasurement*>(measurement_ptr);
    return Json::Value(jsoncpp_real_cast(rmpm->get()));
}

static Json::Value renyi_sign_walk_measurement_json_repr (const Measurement<RenyiSignWalk> *measurement_ptr)
{
    const RenyiSignMeasurement *rsm = boost::polymorphic_downcast<const RenyiSignMeasurement*>(measurement_ptr);
    return complex_to_json_array(rsm->get());
}

template <unsigned int DIM>
static int do_simulation (const Json::Value &json_input, rng_class &rng);

int main (int argc, char *argv[])
{
    // take json input and perform a simulation

    Json::Value json_input;

    if (argc < 2) {
        // read from STDIN
        std::cin >> json_input;
    } else if (argc == 2) {
        // read from given file
        std::ifstream input_file(argv[1]);
        if (!input_file.good()) {
            std::cerr << "cannot open file: " << argv[1] << std::endl;
            return 1;
        }
        input_file >> json_input;
        input_file.close();
    } else {
        std::cerr << "invalid number of commandline arguments" << std::endl;
        return 1;
    }

    try {
        ensure_object(json_input);
        const char * const json_input_required[] = { "rng", "system", "simulation", NULL };
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
        const char * const json_system_required[] = { "lattice", "wavefunction", NULL };
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
    // finish setting up the lattice
    const Json::Value &json_lattice_size = json_input["system"]["lattice"]["size"];
    boost::array<int, DIM> lattice_size_array;
    for (unsigned int i = 0; i < DIM; ++i)
        lattice_size_array[i] = json_lattice_size[i].asInt();
    const boost::shared_ptr<const NDLattice<DIM> > lattice(new NDLattice<DIM>(lattice_size_array));

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
        boost::shared_ptr<const OrbitalDefinitions> orbitals = parse_json_orbitals<DIM>(json_wavefunction["orbitals"], lattice);
        std::vector<std::vector<unsigned int> > filling;
        filling.push_back(some_random_filling<DIM>(orbitals->get_N_filled(), *lattice, rng));
        wf.reset(new FreeFermionWavefunctionAmplitude(PositionArguments(filling, lattice->total_sites()), orbitals));
    } else if (std::strcmp(json_wavefunction_type_cstr, "dbl") == 0) {
        // dbl wavefunction
        const char * const json_dbl_wavefunction_required[] = { "type", "orbitals-d1", "orbitals-d2", NULL };
        const char * const json_dbl_wavefunction_allowed[] = { "type", "orbitals-d1", "orbitals-d2", "exponent-d1", "exponent-d2", NULL };
        ensure_required(json_wavefunction, json_dbl_wavefunction_required);
        ensure_only(json_wavefunction, json_dbl_wavefunction_allowed);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d1 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d1"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d2 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d2"], lattice);
        if (orbitals_d1->get_N_filled() != orbitals_d2->get_N_filled())
            throw ParseError("d1 and d2 have different number of orbitals");
        std::vector<std::vector<unsigned int> > filling;
        filling.push_back(some_random_filling<DIM>(orbitals_d1->get_N_filled(), *lattice, rng));
        wf.reset(new DBLWavefunctionAmplitude(PositionArguments(filling, lattice->total_sites()),
                                              orbitals_d1, orbitals_d2,
                                              json_get_double(json_wavefunction, "exponent-d1", 1.0),
                                              json_get_double(json_wavefunction, "exponent-d2", 1.0)));
    } else if (std::strcmp(json_wavefunction_type_cstr, "dmetal") == 0) {
        // dmetal wavefunction
        const char * const json_dmetal_wavefunction_required[] = { "type", "orbitals-d1", "orbitals-d2", "orbitals-f_up", "orbitals-f_down", NULL };
        const char * const json_dmetal_wavefunction_allowed[] = { "type", "orbitals-d1", "orbitals-d2", "orbitals-f_up", "orbitals-f_down", "exponent-d1", "exponent-d2", "exponent-f_up", "exponent-f_down", NULL };
        ensure_required(json_wavefunction, json_dmetal_wavefunction_required);
        ensure_only(json_wavefunction, json_dmetal_wavefunction_allowed);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d1 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d1"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d2 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d2"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_f_up = parse_json_orbitals<DIM>(json_wavefunction["orbitals-f_up"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_f_down = parse_json_orbitals<DIM>(json_wavefunction["orbitals-f_down"], lattice);
        if (orbitals_d1->get_N_filled() != orbitals_d2->get_N_filled())
            throw ParseError("d1 and d2 have different number of orbitals");
        if (orbitals_f_up->get_N_filled() + orbitals_f_down->get_N_filled() != orbitals_d1->get_N_filled())
            throw ParseError("number of orbitals in f_up and f_down must sum to number of orbitals in d1");
        std::vector<std::vector<unsigned int> > filling;
        filling.push_back(some_random_filling<DIM>(orbitals_f_up->get_N_filled(), *lattice, rng));
        filling.push_back(some_random_filling<DIM>(orbitals_f_down->get_N_filled(), *lattice, rng));
        wf.reset(new DMetalWavefunctionAmplitude(PositionArguments(filling, lattice->total_sites()),
                                                 orbitals_d1, orbitals_d2, orbitals_f_up, orbitals_f_down,
                                                 json_get_double(json_wavefunction, "exponent-d1", 1.0),
                                                 json_get_double(json_wavefunction, "exponent-d2", 1.0),
                                                 json_get_double(json_wavefunction, "exponent-f_up", 1.0),
                                                 json_get_double(json_wavefunction, "exponent-f_down", 1.0)));
    } else {
        throw ParseError("invalid wavefunction type");
    }

    // begin setting up the simulation
    const Json::Value &json_simulation = json_input["simulation"];
    ensure_object(json_simulation);
    const char * const json_simulation_required[] = { "walk-type", "measurement-steps", NULL };
    ensure_required(json_simulation, json_simulation_required);
#define JSON_SIMULATION_GLOBAL_ALLOWED "walk-type", "measurement-steps", "equilibrium-steps", "initial-positions"

    // determine how many steps to do
    unsigned int measurement_steps, equilibrium_steps = 0;
    const Json::Value &json_measurement_steps = json_simulation["measurement-steps"];
    if (!(json_measurement_steps.isIntegral() && json_measurement_steps.asInt() > 0))
        throw ParseError("measurement-steps must be a positive integer");
    measurement_steps = json_measurement_steps.asUInt();
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
    Json::Value json_measurement_output(Json::arrayValue);
    Json::Value json_final_positions_output, json_monte_carlo_stats_output;
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
            bool success = search_for_filling_with_nonzero_amplitude<DIM>(*wf, *lattice, rng);
            if (!success)
                throw ParseError("could not find a filling with nonzero amplitude");
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<Measurement<StandardWalk> > > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i)
            measurements.push_back(parse_standard_walk_measurement_definition<DIM>(json_measurements[i], *wf, *lattice));

        // set up and perform walk
        StandardWalk walk(wf);
        MetropolisSimulation<StandardWalk> sim(walk, measurements, equilibrium_steps, rng());
        sim.iterate(measurement_steps);

        // store json
        for (std::list<boost::shared_ptr<Measurement<StandardWalk> > >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(standard_walk_measurement_json_repr<DIM>(i->get(), walk.get_wavefunction()));
        }
        json_final_positions_output = positions_json_repr(sim.get_walk().get_wavefunction().get_positions());
        json_monte_carlo_stats_output = monte_carlo_stats_json_repr(sim);

    } else if (std::strcmp(json_walk_type_cstr, "renyi-mod") == 0) {

        // RENYI MOD WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "measurements", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "measurements", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up initial positions of wavefunctions
        boost::shared_ptr<WavefunctionAmplitude> wf2(wf->clone());
        if (json_simulation.isMember("initial-positions")) {
            ensure_array(json_simulation["initial-positions"], 2);
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"][0u]);
            set_wavefunction_positions_from_json(*wf2, json_simulation["initial-positions"][1u]);
        } else {
            bool success = (search_for_filling_with_nonzero_amplitude<DIM>(*wf, *lattice, rng)
                            && search_for_filling_with_nonzero_amplitude<DIM>(*wf2, *lattice, rng));
            if (!success)
                throw ParseError("could not find a filling with nonzero amplitude");
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<Measurement<RenyiModWalk> > > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i)
            measurements.push_back(parse_renyi_mod_walk_measurement_definition<DIM>(json_measurements[i], *lattice));

        // set up and perform walk
        RenyiModWalk walk(wf, wf2);
        MetropolisSimulation<RenyiModWalk> sim(walk, measurements, equilibrium_steps, rng());
        sim.iterate(measurement_steps);

        // store json
        for (std::list<boost::shared_ptr<Measurement<RenyiModWalk> > >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(renyi_mod_walk_measurement_json_repr(i->get()));
        }
        json_final_positions_output = Json::Value(Json::arrayValue);
        json_final_positions_output.append(positions_json_repr(sim.get_walk().get_phialpha1().get_positions()));
        json_final_positions_output.append(positions_json_repr(sim.get_walk().get_phialpha2().get_positions()));
        json_monte_carlo_stats_output = monte_carlo_stats_json_repr(sim);

    } else if (std::strcmp(json_walk_type_cstr, "renyi-mod/possible") == 0) {

        // RENYI MOD/POSSIBLE WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", "measurements", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "subsystem", "measurements", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem<DIM>(json_simulation["subsystem"], *lattice));

        // set up initial positions of wavefunctions
        boost::shared_ptr<WavefunctionAmplitude> wf2(wf->clone());
        if (json_simulation.isMember("initial-positions")) {
            ensure_array(json_simulation["initial-positions"], 2);
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"][0u]);
            set_wavefunction_positions_from_json(*wf2, json_simulation["initial-positions"][1u]);
            if (!count_subsystem_particle_counts_for_match(*wf, *wf2, *subsystem))
                throw ParseError("The initial positions of each copy must have the same numbers/types of particles in the subsystem");
        } else {
            bool success = search_for_filling_with_nonzero_amplitude<DIM>(*wf, *lattice, rng);
            if (!success)
                throw ParseError("could not find a filling with nonzero amplitude");
            // We need two copies of the system, each of which has the same
            // number of particles in the subsystem.  So for now we just
            // initialize both copies with the same exact positions.
            wf2 = wf;
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<Measurement<RenyiModPossibleWalk> > > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i)
            measurements.push_back(parse_renyi_mod_possible_walk_measurement_definition<DIM>(json_measurements[i]));

        // set up and perform walk
        RenyiModPossibleWalk walk(wf, wf2, subsystem);
        MetropolisSimulation<RenyiModPossibleWalk> sim(walk, measurements, equilibrium_steps, rng());
        sim.iterate(measurement_steps);

        // store json
        for (std::list<boost::shared_ptr<Measurement<RenyiModPossibleWalk> > >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(renyi_mod_possible_walk_measurement_json_repr(i->get()));
        }
        json_final_positions_output = Json::Value(Json::arrayValue);
        json_final_positions_output.append(positions_json_repr(sim.get_walk().get_phialpha1().get_positions()));
        json_final_positions_output.append(positions_json_repr(sim.get_walk().get_phialpha2().get_positions()));
        json_monte_carlo_stats_output = monte_carlo_stats_json_repr(sim);

    } else if (std::strcmp(json_walk_type_cstr, "renyi-sign") == 0) {

        // RENYI SIGN WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", "measurements", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "subsystem", "measurements", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem<DIM>(json_simulation["subsystem"], *lattice));

        // set up initial positions of wavefunctions
        boost::shared_ptr<WavefunctionAmplitude> wf2(wf->clone());
        if (json_simulation.isMember("initial-positions")) {
            ensure_array(json_simulation["initial-positions"], 2);
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"][0u]);
            set_wavefunction_positions_from_json(*wf2, json_simulation["initial-positions"][1u]);
            if (!count_subsystem_particle_counts_for_match(*wf, *wf2, *subsystem))
                throw ParseError("The initial positions of each copy must have the same numbers/types of particles in the subsystem");
        } else {
            bool success = search_for_filling_with_nonzero_amplitude<DIM>(*wf, *lattice, rng);
            if (!success)
                throw ParseError("could not find a filling with nonzero amplitude");
            // We need two copies of the system, each of which has the same
            // number of particles in the subsystem.  So for now we just
            // initialize both copies with the same exact positions.
            wf2 = wf;
        }

        // set up measurement(s)
        const Json::Value &json_measurements = json_simulation["measurements"];
        ensure_array(json_measurements);
        std::list<boost::shared_ptr<Measurement<RenyiSignWalk> > > measurements;
        for (unsigned int i = 0; i < json_measurements.size(); ++i)
            measurements.push_back(parse_renyi_sign_walk_measurement_definition<DIM>(json_measurements[i]));

        // set up and perform walk
        RenyiSignWalk walk(wf, wf2, subsystem);
        MetropolisSimulation<RenyiSignWalk> sim(walk, measurements, equilibrium_steps, rng());
        sim.iterate(measurement_steps);

        // store json
        for (std::list<boost::shared_ptr<Measurement<RenyiSignWalk> > >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(renyi_sign_walk_measurement_json_repr(i->get()));
        }
        json_final_positions_output = Json::Value(Json::arrayValue);
        json_final_positions_output.append(positions_json_repr(sim.get_walk().get_phialpha1().get_positions()));
        json_final_positions_output.append(positions_json_repr(sim.get_walk().get_phialpha2().get_positions()));
        json_monte_carlo_stats_output = monte_carlo_stats_json_repr(sim);

    } else {
        throw ParseError("invalid walk type");
    }

    Json::Value json_output(Json::objectValue);
    json_output["final-positions"] = json_final_positions_output;
    json_output["measurements"] = json_measurement_output;
    json_output["run-information"] = run_information.json_info();
    json_output["monte-carlo-stats"] = json_monte_carlo_stats_output;
    std::cout << json_output;

    return 0;
}
