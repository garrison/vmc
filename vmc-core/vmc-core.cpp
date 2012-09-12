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

boost::shared_ptr<Wavefunction::Amplitude> create_wfa_from_json (const char *json_input_str, const boost::shared_ptr<const Lattice> &lattice, std::auto_ptr<RandomNumberGenerator> &rng)
{
    Json::Value json_system;
    {
        std::istringstream iss(json_input_str, std::istringstream::in);
        iss >> json_system;
    }

    // set up the wavefunction
    ensure_object(json_system);
    const char * const json_system_required[] = { "wavefunction", NULL };
    ensure_required(json_system, json_system_required);
    ensure_only(json_system, json_system_required);
    boost::shared_ptr<Wavefunction> wf(create_wavefunction(json_system["wavefunction"], lattice));

    // set up the Wavefunction::Amplitude
    boost::shared_ptr<Wavefunction::Amplitude> wfa(wf->create_nonzero_wavefunctionamplitude(wf, *rng));
    if (!wfa)
        throw ParseError("could not find a nonzero wavefunction amplitude");

    return wfa;
}
