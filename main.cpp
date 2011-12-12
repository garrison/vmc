#include <iostream>
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

#include "NDLattice.hpp"
#include "random-combination.hpp"
#include "PositionArguments.hpp"
#include "FilledOrbitals.hpp"
#include "FreeFermionWavefunctionAmplitude.hpp"
#include "DBLWavefunctionAmplitude.hpp"
#include "DBMWavefunctionAmplitude.hpp"
#include "MetropolisSimulation.hpp"
#include "StandardWalk.hpp"
#include "DensityDensityMeasurement.hpp"
#include "GreenMeasurement.hpp"
#include "RenyiModWalk.hpp"
#include "RenyiModMeasurement.hpp"
#include "RenyiSignWalk.hpp"
#include "RenyiSignMeasurement.hpp"
#include "SimpleSubsystem.hpp"

// fixme: there is currently no reason for this to be a template
// fixme: this function is more general than main.cpp
template <unsigned int DIM>
static PositionArguments some_random_filling (unsigned int N_filled, const NDLattice<DIM> &lattice, rng_class &rng)
{
    std::vector<unsigned int> v;
    random_combination(v, N_filled, lattice.total_sites(), rng);
    return PositionArguments(v, lattice.total_sites());
}

// fixme: there is currently no reason for this to be a template
// fixme: this function is more general than main.cpp
template <unsigned int DIM>
static bool search_for_filling_with_nonzero_amplitude (WavefunctionAmplitude &wf, const NDLattice<DIM> &lattice, rng_class &rng)
{
    unsigned int attempts = 1; // assume that one attempt has already been completed
    while (wf.psi() == amplitude_t(0)) {
        if (attempts++ == 10000)
            return false;
        wf.reset(some_random_filling<DIM>(wf.get_positions().get_N_filled(), lattice, rng));
    }
    return true;
}

// fixme: this function is more general than main.cpp
static unsigned int count_N_subsystem (const WavefunctionAmplitude &wf, const Subsystem &subsystem)
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

template <unsigned int DIM>
boost::shared_ptr<const OrbitalDefinitions> parse_json_orbitals (const Json::Value &json_orbitals, const boost::shared_ptr<const NDLattice<DIM> > &lattice)
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
boost::shared_ptr<Measurement<StandardWalk> > parse_standard_walk_measurement_definition (const Json::Value &json_measurement_def)
{
    ensure_object_with_type_field_as_string(json_measurement_def);
    if (std::strcmp(json_measurement_def["type"].asCString(), "density-density") == 0) {
        const char * const json_density_density_allowed[] = { "type", NULL };
        ensure_only(json_measurement_def, json_density_density_allowed);
        return boost::make_shared<DensityDensityMeasurement<DIM> >();
    } else if (std::strcmp(json_measurement_def["type"].asCString(), "green") == 0) {
        const char * const json_density_density_allowed[] = { "type", NULL };
        ensure_only(json_measurement_def, json_density_density_allowed);
        return boost::make_shared<GreenMeasurement<DIM> >();
    } else {
        throw ParseError("invalid standard walk measurement type");
    }
}

template <unsigned int DIM>
boost::shared_ptr<Measurement<RenyiModWalk> > parse_renyi_mod_walk_measurement_definition (const Json::Value &json_measurement_def, const NDLattice<DIM> &lattice)
{
    boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem<DIM>(json_measurement_def, lattice));
    return boost::make_shared<RenyiModMeasurement>(subsystem);
}

static void set_wavefunction_positions_from_json (WavefunctionAmplitude &wf, const Json::Value &json_initial_positions)
{
    const unsigned int N_filled = wf.get_positions().get_N_filled();
    const unsigned int N_sites = wf.get_positions().get_N_sites();
    ensure_array(json_initial_positions, N_filled);
    std::vector<unsigned int> v(N_filled);
    std::set<unsigned int> vs;
    for (unsigned int i = 0; i < N_filled; ++i) {
        const Json::Value &json_pos = json_initial_positions[i];
        if (!json_pos.isIntegral())
            throw ParseError("expecting integer");
        if (!(json_pos.asInt() >= 0 && json_pos.asUInt() < N_sites))
            throw ParseError("invalid site index");
        bool inserted = vs.insert(json_pos.asUInt()).second;
        if (!inserted)
            throw ParseError("position index specified twice, but double occupancy is not allowed");
        v[i] = json_pos.asUInt();
    }
    wf.reset(PositionArguments(v, N_sites));
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
static Json::Value standard_walk_measurement_json_repr (const Measurement<StandardWalk> *measurement_ptr)
{
    const DensityDensityMeasurement<DIM> *ddm = dynamic_cast<const DensityDensityMeasurement<DIM>*>(measurement_ptr);
    const GreenMeasurement<DIM> *gm = dynamic_cast<const GreenMeasurement<DIM>*>(measurement_ptr);
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

static Json::Value renyi_sign_walk_measurement_json_repr (const RenyiSignMeasurement *measurement)
{
    return complex_to_json_array(measurement->get());
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
        wf.reset(new FreeFermionWavefunctionAmplitude(some_random_filling<DIM>(orbitals->get_N_filled(), *lattice, rng), orbitals));
    } else if (std::strcmp(json_wavefunction_type_cstr, "dbl") == 0) {
        // dbl wavefunction
        const char * const json_dbl_wavefunction_required[] = { "type", "orbitals-d1", "orbitals-d2", NULL };
        ensure_required(json_wavefunction, json_dbl_wavefunction_required);
        ensure_only(json_wavefunction, json_dbl_wavefunction_required);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d1 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d1"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d2 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d2"], lattice);
        if (orbitals_d1->get_N_filled() != orbitals_d2->get_N_filled())
            throw ParseError("d1 and d2 have different number of orbitals");
        wf.reset(new DBLWavefunctionAmplitude(some_random_filling<DIM>(orbitals_d1->get_N_filled(), *lattice, rng),
                                              orbitals_d1, orbitals_d2));
    } else if (std::strcmp(json_wavefunction_type_cstr, "dbm") == 0) {
        // dbm wavefunction
        const char * const json_dbm_wavefunction_required[] = { "type", "orbitals-d1", "orbitals-d2", "orbitals-f_up", "orbitals-f_down", NULL };
        ensure_required(json_wavefunction, json_dbm_wavefunction_required);
        ensure_only(json_wavefunction, json_dbm_wavefunction_required);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d1 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d1"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_d2 = parse_json_orbitals<DIM>(json_wavefunction["orbitals-d2"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_f_up = parse_json_orbitals<DIM>(json_wavefunction["orbitals-f_up"], lattice);
        boost::shared_ptr<const OrbitalDefinitions> orbitals_f_down = parse_json_orbitals<DIM>(json_wavefunction["orbitals-f_down"], lattice);
        if (orbitals_d1->get_N_filled() != orbitals_d2->get_N_filled())
            throw ParseError("d1 and d2 have different number of orbitals");
        if (orbitals_f_up->get_N_filled() + orbitals_f_down->get_N_filled() != orbitals_d1->get_N_filled())
            throw ParseError("number of orbitals in f_up and f_down must sum to number of orbitals in d1");
        wf.reset(new DBMWavefunctionAmplitude(some_random_filling<DIM>(orbitals_d1->get_N_filled(), *lattice, rng),
                                              orbitals_d1, orbitals_d2, orbitals_f_up, orbitals_f_down));
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
            measurements.push_back(parse_standard_walk_measurement_definition<DIM>(json_measurements[i]));

        // set up and perform walk
        StandardWalk walk(wf);
        MetropolisSimulation<StandardWalk> sim(walk, measurements, equilibrium_steps, rng());
        sim.iterate(measurement_steps);

        // store json
        std::cerr << "accepted " << (100.0 * sim.steps_accepted() / sim.steps_completed()) << "%" << std::endl;
        // fixme: add current positions to json as well as steps accepted/rejected/totally-rejected, etc
        for (std::list<boost::shared_ptr<Measurement<StandardWalk> > >::const_iterator i = measurements.begin(); i != measurements.end(); ++i) {
            json_measurement_output.append(standard_walk_measurement_json_repr<DIM>(i->get()));
        }

    } else if (std::strcmp(json_walk_type_cstr, "renyi-mod") == 0) {

        // RENYI MOD WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "measurement-subsystems", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "measurement-subsystems", NULL };
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
        const Json::Value &json_measurements = json_simulation["measurement-subsystems"];
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

    } else if (std::strcmp(json_walk_type_cstr, "renyi-sign") == 0) {

        // RENYI SIGN WALK

        // ensure correct json properties are given
        const char * const json_simulation_additional_required[] = { "subsystem", NULL };
        ensure_required(json_simulation, json_simulation_additional_required);
        const char * const json_simulation_only[] = { JSON_SIMULATION_GLOBAL_ALLOWED, "subsystem", NULL };
        ensure_only(json_simulation, json_simulation_only);

        // set up subsystem
        boost::shared_ptr<const Subsystem> subsystem(parse_json_subsystem<DIM>(json_simulation["subsystem"], *lattice));

        // set up initial positions of wavefunctions
        boost::shared_ptr<WavefunctionAmplitude> wf2(wf->clone());
        if (json_simulation.isMember("initial-positions")) {
            ensure_array(json_simulation["initial-positions"], 2);
            set_wavefunction_positions_from_json(*wf, json_simulation["initial-positions"][0u]);
            set_wavefunction_positions_from_json(*wf2, json_simulation["initial-positions"][1u]);
            if (count_N_subsystem(*wf, *subsystem) != count_N_subsystem(*wf2, *subsystem))
                throw ParseError("The initial positions of each copy must have the same number of particles in the subsystem");
        } else {
            bool success = search_for_filling_with_nonzero_amplitude<DIM>(*wf, *lattice, rng);
            if (!success)
                throw ParseError("could not find a filling with nonzero amplitude");
            // We need two copies of the system, each of which has the same
            // number of particles in the subsystem.  So for now we just
            // initialize both copies with the same exact positions.
            wf2 = wf;
        }

        // set up measurement
        boost::shared_ptr<RenyiSignMeasurement> measurement(new RenyiSignMeasurement);

        // set up and perform walk
        RenyiSignWalk walk(wf, wf2, subsystem);
        MetropolisSimulation<RenyiSignWalk> sim(walk, measurement, equilibrium_steps, rng());
        sim.iterate(measurement_steps);

        // store json
        json_measurement_output.append(renyi_sign_walk_measurement_json_repr(measurement.get()));

    } else {
        throw ParseError("invalid walk type");
    }

    std::cout << json_measurement_output;

    return 0;
}
