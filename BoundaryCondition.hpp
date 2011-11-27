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
    // the argument specifies what fraction (of $2\pi$) the phase is increased
    // when moving once through the system in a given direction.  1 corresponds
    // to periodic boundary conditions; 2 corresponds to antiperiodic; etc.
    explicit BoundaryCondition (const boost::rational<int> &p_)
        : m_p(p_),
          m_p_floor(p_ == 1 ? 0 : m_p),
          m_phase(calculate_phase(p_))
        {
            BOOST_ASSERT(p_ > 0 && p_ <= 1);
        }

    // the argument specifies how many times one must move through the system
    // to return to the same phase.  1 means periodic; 2 means antiperiodic;
    // etc.  this constructor is given as a convenience, as it allows us to
    // specify boundary conditions with a single integer in any case where they
    // are not crazy.
    explicit BoundaryCondition (unsigned int p_)
        : m_p(boost::rational<int>(1, p_)),
          m_p_floor(p_ == 1 ? 0 : m_p),
          m_phase(calculate_phase(p_))
        {
        }

    BoundaryCondition (void) // uninitialized default constructor
        : m_p(0),
          m_p_floor(0),
          m_phase(0)
        {
        }

    // returns a value in (0, 1]
    boost::rational<int> p (void) const
        {
            BOOST_ASSERT(m_p != 0); // otherwise it is uninitialized
            return m_p;
        }

    // returns a value in [0, 1)
    boost::rational<int> p_floor (void) const
        {
            BOOST_ASSERT(m_p != 0); // otherwise it is uninitialized
            return m_p_floor;
        }

    // returns the phase change when one crosses the boundary in the positive
    // direction
    phase_t phase (void) const
        {
            BOOST_ASSERT(m_p != 0); // otherwise it is uninitialized
            return m_phase;
        }

private:
    static phase_t calculate_phase (const boost::rational<int> &p)
        {
            // if we can return an exact value, do so
            if (p == boost::rational<int>(1))
                return phase_t(1, 0);
            else if (p == boost::rational<int>(1, 2))
                return phase_t(-1, 0);
            else if (p == boost::rational<int>(1, 4))
                return phase_t(0, 1);
            else if (p == boost::rational<int>(3, 4))
                return phase_t(0, -1);
            else
                // if not, fall back using the exponential function
                return std::exp(complex_t(0, 1) * complex_t(2 * boost::math::constants::pi<real_t>() * boost::rational_cast<real_t>(p)));
        }

    boost::rational<int> m_p, m_p_floor;
    phase_t m_phase;
};

// fixme: initialize these once only?
static const BoundaryCondition periodic_bc(1);
static const BoundaryCondition antiperiodic_bc(2);

#endif
