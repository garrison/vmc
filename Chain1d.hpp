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

    Chain1d::amplitude_t phi(int n, int r) const
	{
	    BOOST_ASSERT(n >= 0 && n < N_filled);
	    BOOST_ASSERT(r >= 0 && r < N_sites);

	    const double two_pi = 2 * boost::math::constants::pi<double>();
	    int kbar = (n / 2 + 1) * ((n & 1) * -2 + 1); // fill each k value in order, alternating +k, -k
	    const std::complex<double> im_unit(0, 1);
	    return exp(im_unit * std::complex<double>(two_pi * kbar * r / N_sites));
	}

private:
    int N_filled, N_sites;

private:
    // disable default constructor
    Chain1d(void);
};

class Chain1dArguments
// position arguments, and knowledge of constraints on them, are contained in this class.
{
private:
    //boost::shared_ptr<Chain1d> wf;
    std::vector<Chain1d::position_t> r;
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
	    r[v] = position;
	}

    void swap_positions (int v1, int v2)
	{
	    BOOST_ASSERT(v1 >= 0 && v1 < (int)r.size());
	    BOOST_ASSERT(v2 >= 0 && v2 < (int)r.size());
	    std::swap(r[v1], r[v2]);
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

class Chain1dRenyiWalk
{
public:
    enum RenyiWalkType {
	SWAPA_MOD,
	SWAPA_SIGN
    };
private:
    boost::shared_ptr<Chain1d> wf;
    Chain1dArguments r1, r2;
    CeperlyMatrix<Chain1d::amplitude_t> phialpha1, phialpha2, phibeta1, phibeta2;
    const Subsystem<Chain1d> *subsystem;
    int N_subsystem1, N_subsystem2;
    RenyiWalkType walk_type;
    bool transition_in_progress;
public:
    Chain1dRenyiWalk (const Chain1d &wf_, const Subsystem<Chain1d> *subsystem_, RenyiWalkType walk_type_, rng_class &rng);
    //Chain1dRenyiWalk (const Chain1d &wf_, const Subsystem<Chain1d> *subsystem_, RenyiWalkType walk_type_, const std::vector<Chain1d::position_t> &r);
    probability_t compute_probability_ratio_of_random_transition (rng_class &rng);
    void accept_transition (void);
private:
    // disable the default constructor
    Chain1dRenyiWalk (void);
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

#endif
