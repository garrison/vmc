#ifndef _VMC_BOUNDARY_CONDITION_HPP
#define _VMC_BOUNDARY_CONDITION_HPP

#include <cmath>

#include <boost/assert.hpp>
#include <boost/rational.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_complex.hpp>

#include "lw_vector.hpp"
#include "vmc-typedefs.hpp"

/**
 * Represents a boundary condition in one dimension for a system that exists on
 * an N-dimensional torus.  Both periodic and antiperiodic boundary conditions
 * are supported, as well as a variety of "twisted" boundary conditions, in
 * which the relevant complex quantity advances by some arbitrary phase (given
 * by a rational number) when one wraps around the system a single time.
 */
template <typename PhaseType>
class BoundaryCondition
{
public:
    /**
     * Constructor.
     *
     * @param p_ specifies what fraction (of \f$2\pi\f$) the phase is increased
     * when moving once through the system in the relevant direction.  1
     * corresponds to periodic boundary conditions; 1/2 corresponds to
     * antiperiodic; etc.  0 corresponds to open boundary conditions.
     */
    explicit BoundaryCondition (const boost::rational<int> &p_)
        : m_p(p_),
          m_phase(calculate_phase<PhaseType>(p_))
        {
            BOOST_ASSERT(p_ >= 0 && p_ <= 1);
        }

    /**
     * Uninitialized default constructor
     */
    BoundaryCondition (void)
        : m_p(-1),
          m_phase(0)
        {
        }

    /**
     * Returns a value in [0, 1]
     */
    boost::rational<int> p (void) const
        {
            BOOST_ASSERT(m_p != -1); // otherwise it is uninitialized
            return m_p;
        }

    /**
     * Returns the phase change when one crosses the boundary in the positive
     * direction.  This will be zero for open boundary conditions, or will be
     * along the unit circle for any type of periodic boundary conditions.
     */
    PhaseType phase (void) const
        {
            BOOST_ASSERT(m_p != -1); // otherwise it is uninitialized
            return m_phase;
        }

    bool is_initialized (void) const
        {
            return m_p != -1;
        }

    bool operator== (const BoundaryCondition &other) const
        {
            return (m_p == other.m_p);
        }

    bool operator!= (const BoundaryCondition &other) const
        {
            return (m_p != other.m_p);
        }

private:
    /**
     * This function is called to initialize the data member m_phase during
     * object construction.  We are unsure whether PhaseType is real or complex,
     * so we implement both and the compiler chooses the correct one using
     * template specialization.
     */
    template <typename PHASE_T>
    static typename boost::enable_if<boost::is_complex<PHASE_T>, PhaseType>::type calculate_phase (const boost::rational<int> &p)
        {
            // if we can return an exact value, do so
            if (p == 0)
                return PhaseType(0); // open
            else if (p == boost::rational<int>(1))
                return PhaseType(1); // periodic
            else if (p == boost::rational<int>(1, 2))
                return PhaseType(-1); // antiperiodic
            else if (p == boost::rational<int>(1, 4))
                return PhaseType(0, 1);
            else if (p == boost::rational<int>(3, 4))
                return PhaseType(0, -1);
            else
                // cannot return an exact value, so fall back using the exponential function
                return std::exp(complex_t(0, 1) * complex_t(2 * boost::math::constants::pi<real_t>() * boost::rational_cast<real_t>(p)));
        }

    template <typename PHASE_T>
    static typename boost::disable_if<boost::is_complex<PHASE_T>, PhaseType>::type calculate_phase (const boost::rational<int> &p)
        {
            if (p == boost::rational<int>(1))
                return PhaseType(1); // periodic
            else if (p == boost::rational<int>(1, 2))
                return PhaseType(-1); // antiperiodic

            // the only remaining possibility for a real PhaseType is open boundary
            // conditions
            BOOST_ASSERT(p == 0);
            return PhaseType(0);
        }

    boost::rational<int> m_p;
    PhaseType m_phase;

public:
    static const BoundaryCondition<PhaseType> open;
    static const BoundaryCondition<PhaseType> periodic;
    static const BoundaryCondition<PhaseType> antiperiodic;
};

/** boundary conditions in each direction */
template <typename PhaseType>
using BoundaryConditions = lw_vector<BoundaryCondition<PhaseType>, MAX_DIMENSION>;

// named boundary conditions
template <typename PhaseType>
const BoundaryCondition<PhaseType> BoundaryCondition<PhaseType>::open = BoundaryCondition<PhaseType>(boost::rational<int>(0));

template <typename PhaseType>
const BoundaryCondition<PhaseType> BoundaryCondition<PhaseType>::periodic = BoundaryCondition<PhaseType>(boost::rational<int>(1));

template <typename PhaseType>
const BoundaryCondition<PhaseType> BoundaryCondition<PhaseType>::antiperiodic = BoundaryCondition<PhaseType>(boost::rational<int>(1, 2));

#endif
