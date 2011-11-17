#ifndef _BOUNDARY_CONDITION_HPP
#define _BOUNDARY_CONDITION_HPP

#include <cmath>
#include <cstddef>

#include <boost/assert.hpp>
#include <boost/rational.hpp>
#include <boost/math/constants/constants.hpp>

#include "vmc-typedefs.hpp"

// fixme: make it not be inline?
class BoundaryCondition
{
public:
    explicit BoundaryCondition (const boost::rational<unsigned int> &p_)
	: m_p(p_),
	  m_phase(calculate_phase(p_))
	{
	    BOOST_ASSERT(p_ > 0 && p_ <= 1);
	}

    explicit BoundaryCondition (unsigned int p_)
	: m_p(boost::rational<unsigned int>(1, p_)),
	  m_phase(calculate_phase(m_p))
	{
	}

    boost::rational<unsigned int> p (void) const
	{
	    return m_p;
	}

    phase_t phase (void) const
	{
	    return m_phase;
	}

private:
    static phase_t calculate_phase (const boost::rational<unsigned int> &p)
	{
	    // if we can return an exact value, do so
	    if (p == boost::rational<unsigned int>(1))
		return phase_t(1, 0);
	    else if (p == boost::rational<unsigned int>(1, 2))
		return phase_t(-1, 0);
	    else if (p == boost::rational<unsigned int>(1, 4))
		return phase_t(0, 1);
	    else if (p == boost::rational<unsigned int>(3, 4))
		return phase_t(0, -1);
	    else
		// if not, fall back using the exponential function
		return std::exp(complex_t(0, 1) * complex_t(2 * boost::math::constants::pi<real_t>() * boost::rational_cast<real_t>(p)));
	}

    boost::rational<unsigned int> m_p;
    phase_t m_phase;
};

// fixme: initialize these once only?
static const BoundaryCondition periodic_bc(1);
static const BoundaryCondition antiperiodic_bc(2);

#endif
