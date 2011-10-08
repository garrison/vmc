#ifndef _CHAIN1D_HPP
#define _CHAIN1D_HPP

#include <vector>
#include <memory>
#include <cmath>

#include <boost/assert.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/math/constants/constants.hpp>

#include "vmc-typedefs.hpp"
#include "CeperlyMatrix.hpp"

// contains all functions and parameters of the wave function itself.
// Generally the actual computation of \Psi(R) will occur within the Walk
// class, since things can be computed less redundantly that way.
class Chain1d
{
public:
    typedef int position_t;
    typedef complex_amplitude_t amplitude_t;

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

	    const double pi = boost::math::constants::pi<double>();
	    const std::complex<double> im_unit(0, 1);
	    return exp(im_unit * std::complex<double>(2 * pi * r * (n + 1) / N_sites));
	}

private:
    int N_filled, N_sites;

private:
    // disable default constructor
    Chain1d(void);
};

class Chain1dWalk // Both the arguments themselves (and knowledge of
		       // constraints on them) as well as the methods for
		       // performing the random walk are in this class.  In
		       // theory they should be in two different classes, but
		       // we wouldn't gain anything from that (and the concepts
		       // are closely coupled anyway).
{
private:
    boost::shared_ptr<Chain1d> wf;
    std::vector<Chain1d::position_t> r;
    CeperlyMatrix<Chain1d::amplitude_t> cmat;
    bool transition_in_progress;
public:
    Chain1dWalk (const Chain1d &wf_);
    Chain1dWalk (const Chain1d &wf_, const std::vector<Chain1d::position_t> &r);
    probability_t compute_probability_ratio_of_random_transition (void);
    void finalize_transition (void);
private:
    void initialize_cmat (void);
    // disable the default constructor
    Chain1dWalk (void);
};

#endif
