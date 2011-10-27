#ifndef _CHAIN1D_HPP
#define _CHAIN1D_HPP

#include <algorithm>
#include <vector>
#include <memory>
#include <cmath>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/math/constants/constants.hpp>

#include "vmc-typedefs.hpp"
#include "CeperlyMatrix.hpp"
#include "Subsystem.hpp"
#include "SwappedSystem.hpp"

// NOTE: all "random transitions" are expected to satisfy balance

class Chain1dArguments;

// contains all functions and parameters of the wave function itself.
// Generally the actual computation of \Psi(R) will occur within the Walk
// class, since things can be computed less redundantly that way.
class Chain1d
{
public:
    typedef int position_t;
    typedef complex_amplitude_t amplitude_t;
    typedef Chain1dArguments Arguments;

    Chain1d (int _N_filled, int _N_sites)
	: N_filled(_N_filled),
	  N_sites(_N_sites)
	{
	    BOOST_ASSERT(N_sites > 0 && N_filled > 0);
	    BOOST_ASSERT(N_filled <= N_sites);
	}

    int get_N_sites (void) const
	{
	    return N_sites;
	}

    int get_N_filled (void) const
	{
	    return N_filled;
	}

    Chain1d::amplitude_t phi (int n, int r) const
	{
	    BOOST_ASSERT(n >= 0 && n < N_filled);
	    BOOST_ASSERT(r >= 0 && r < N_sites);

	    const double two_pi = 2 * boost::math::constants::pi<double>();
	    int kbar = ((n + 1) / 2) * ((n & 1) * -2 + 1); // fill each k value in order, alternating +k, -k
	    const std::complex<double> im_unit(0, 1);
	    return exp(im_unit * std::complex<double>(two_pi * kbar * r / N_sites));
	}

private:
    int N_filled, N_sites;

private:
    // disable default constructor
    Chain1d (void);
};

class Chain1dArguments
// position arguments, and knowledge of constraints on them, are contained in this class.
{
private:
    std::vector<Chain1d::position_t> r;
    std::vector<int> positions;
public:
    Chain1dArguments (const Chain1d &wf_, rng_class &rng);
    Chain1dArguments (const Chain1d &wf_, const std::vector<Chain1d::position_t> &r_);

    Chain1d::position_t operator[] (int v) const
	{
	    BOOST_ASSERT(v >= 0 && v < (int)r.size());
	    return r[v];
	}

    size_t size (void) const
	{
	    return r.size();
	}

    void update_position (int v, Chain1d::position_t position)
	{
	    BOOST_ASSERT(v >= 0 && v < (int)r.size());
	    //BOOST_ASSERT(position is valid) // fixme
	    --positions[r[v]];
	    ++positions[position];
	    r[v] = position;
	}

    bool is_occupied (Chain1d::position_t position)
	{
	    return positions[position] != 0;
	}

private:
    // disable the default constructor
    Chain1dArguments (void);
};

class Chain1dWalk
{
private:
    boost::shared_ptr<Chain1d> wf;
    Chain1dArguments r;
    CeperlyMatrix<Chain1d::amplitude_t> cmat;
    bool transition_in_progress;
public:
    Chain1dWalk (const Chain1d &wf_, const Chain1dArguments &arguments_);
    probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
    void accept_transition (void);
private:
    // disable the default constructor
    Chain1dWalk (void);
};

class Chain1dRenyiModWalk
{
private:
    boost::shared_ptr<Chain1d> wf;
    Chain1dArguments r1, r2;
    CeperlyMatrix<Chain1d::amplitude_t> phialpha1, phialpha2;
    std::vector<SwappedSystem<Chain1d> > swapped_system;
    unsigned int transition_copy_in_progress;
    int chosen_particle;
public:
    Chain1dRenyiModWalk (const Chain1d &wf_, const std::vector<Subsystem<Chain1d> *> &subsystems, rng_class &rng);
    probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
    void accept_transition (void);

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phialpha1 (void) const
	{
	    return phialpha1;
	}

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phialpha2 (void) const
	{
	    return phialpha2;
	}

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phibeta1 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index].get_phibeta1();
	}

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phibeta2 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index].get_phibeta2();
	}

    unsigned int get_N_subsystem1 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index].get_N_subsystem1();
	}

    unsigned int get_N_subsystem2 (unsigned int subsystem_index) const
	{
	    BOOST_ASSERT(subsystem_index < swapped_system.size());
	    return swapped_system[subsystem_index].get_N_subsystem2();
	}

    unsigned int get_subsystem_array_size (void) const
	{
	    return swapped_system.size();
	}

private:
    // disable the default constructor
    Chain1dRenyiModWalk (void);
};

class Chain1dRenyiSignWalk
{
private:
    boost::shared_ptr<Chain1d> wf;
    Chain1dArguments r1, r2;
    CeperlyMatrix<Chain1d::amplitude_t> phialpha1, phialpha2;
    SwappedSystem<Chain1d> swapped_system;
    bool transition_in_progress;
public:
    Chain1dRenyiSignWalk (const Chain1d &wf_, const Subsystem<Chain1d> *subsystem_, rng_class &rng);
    probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
    void accept_transition (void);

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phialpha1 (void) const
	{
	    return phialpha1;
	}

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phialpha2 (void) const
	{
	    return phialpha2;
	}

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phibeta1 (void) const
	{
	    return swapped_system.get_phibeta1();
	}

    const CeperlyMatrix<Chain1d::amplitude_t> & get_phibeta2 (void) const
	{
	    return swapped_system.get_phibeta2();
	}

    unsigned int get_N_subsystem1 (void) const
	{
	    return swapped_system.get_N_subsystem1();
	}

    unsigned int get_N_subsystem2 (void) const
	{
	    return swapped_system.get_N_subsystem2();
	}

private:
    // disable the default constructor
    Chain1dRenyiSignWalk (void);
};

class Chain1dContiguousSubsystem : public Subsystem<Chain1d>
{
public:
    Chain1dContiguousSubsystem (Chain1d::position_t subsystem_length_)
	: subsystem_length(subsystem_length_)
	{
	}

    virtual bool particle_is_within (const Chain1d::position_t &position) const
	{
	    BOOST_ASSERT(position >= 0);
	    return position < subsystem_length;
	}

private:
    int subsystem_length;
    // disable default constructor
    Chain1dContiguousSubsystem (void);
};

class Chain1dRenyiModMeasurement
{
public:
    typedef std::vector<double> measurement_value_t;

    Chain1dRenyiModMeasurement (const Chain1dRenyiModWalk &walk)
	: accum(walk.get_subsystem_array_size())
	{
	}

    void measure (const Chain1dRenyiModWalk &walk)
	{
	    for (unsigned int i = 0; i < accum.size(); ++i) {
		if (walk.get_N_subsystem1(i) == walk.get_N_subsystem2(i)) {
		    accum[i] += std::abs(walk.get_phibeta1(i).get_determinant()
					 / walk.get_phialpha1().get_determinant()
					 * walk.get_phibeta2(i).get_determinant()
					 / walk.get_phialpha2().get_determinant());
		}
	    }
	}

    measurement_value_t get (unsigned int measurements_completed) const
	{
	    std::vector<double> rv(accum.size());
	    for (unsigned int i = 0; i < accum.size(); ++i)
		rv[i] = static_cast<double>(accum[i]) / measurements_completed;
	    return rv;
	}

private:
    std::vector<accumulator_t> accum;

    // disable default constructor
    Chain1dRenyiModMeasurement (void);
};

class Chain1dRenyiSignMeasurement
{
    // fixme: in the case of a real wave function, we know that the sign will
    // always evaluate to 1 or -1.  In this case, it may make more sense to
    // store as integers how many measurements were of each value, instead of
    // using a complex<double> as an accumulator.

public:
    typedef Chain1d::amplitude_t measurement_value_t;

    Chain1dRenyiSignMeasurement (const Chain1dRenyiSignWalk &walk)
	: accum(0)
	{
	    (void) walk; // silence warning about unused variable
	}

    void measure (const Chain1dRenyiSignWalk &walk)
	{
	    // we take the argument of each determinant separately instead of
	    // multiplying the determinants together first.  this is necessary
	    // because the determinants tend to be quite large, and multiplying
	    // them can lead to overflow
	    double a = 0;
	    a -= std::arg(walk.get_phialpha1().get_determinant());
	    a -= std::arg(walk.get_phialpha2().get_determinant());
	    a += std::arg(walk.get_phibeta1().get_determinant());
	    a += std::arg(walk.get_phibeta2().get_determinant());
	    const std::complex<double> i(0, 1);
#if 0
	    std::cerr << std::real(std::exp(i * a)) << "   " << std::arg(walk.get_phialpha1().get_determinant())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phialpha2().get_determinant())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phibeta1().get_determinant())/boost::math::constants::pi<double>() << ' ' << std::arg(walk.get_phibeta2().get_determinant())/boost::math::constants::pi<double>() << ' ' << std::endl;
#endif
	    accum += std::exp(i * a);
	}

    measurement_value_t get (unsigned int measurements_completed) const
	{
	    return accum / (real_amplitude_t) measurements_completed;
	}

private:
    measurement_value_t accum; // fixme: complex accumulator_t needed

    // disable default constructor
    Chain1dRenyiSignMeasurement (void);
};

#endif
