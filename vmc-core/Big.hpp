#ifndef _VMC_BIG_HPP
#define _VMC_BIG_HPP

#include <cmath>

#include <boost/assert.hpp>
#include <boost/type_traits/is_complex.hpp>

/**
 * A very big (or small) number, represented as z = A exp(B), where A and B are
 * called the "base" and "exponent", respectively.  B must be a real number,
 * but A can be complex (depending on what T is).
 *
 * This is particularly suited to storing the values of determinants.  The
 * ratios of two Big numbers should be accurate.  Note that this class does not
 * support addition and subtraction, as such operations would be very
 * inaccurate.
 */
template <typename T>
class Big
{
private:
    // The point of the following section of code is to use C++ "partial
    // template specialization" to determine the correct type for the exponent
    template <bool, class>
    class ExponentType;

    // complex type - > real exponent
    template <class T1>
    class ExponentType<true, T1>
    {
    public:
        typedef typename T1::value_type type;
    };

    // real type -> real exponent
    template <class T1>
    class ExponentType<false, T1>
    {
    public:
        typedef T1 type;
    };

public:
    // define result_t, making use of the template specialization above
    typedef typename ExponentType<boost::is_complex<T>::value, T>::type exponent_t;

    /**
     * Explicit creation from a base and exponent
     */
    Big (const T &base, const exponent_t &exponent)
        : m_base(base),
          m_exponent(exponent)
        {
            BOOST_ASSERT(exponent == exponent);
        }

    /*
     * Explicit creation from a T value
     */
    Big (const T &base)
        : m_base(base),
          m_exponent(0)
        {
        }

    /**
     * default constructor, sets it to zero
     */
    Big (void)
        : m_base(0),
          m_exponent(0)
        {
        }

    Big& operator*= (const T &other)
        {
            this->m_base *= other;
            return *this;
        }

    Big<T> operator* (const T &other) const
        {
            return Big<T>(m_base * other, m_exponent);
        }

    Big& operator*= (const Big<T> &other)
        {
            this->m_base *= other.m_base;
            this->m_exponent += other.m_exponent;
            return *this;
        }

    Big<T> operator* (const Big<T> &other) const
        {
            return Big<T>(m_base * other.m_base, m_exponent + other.m_exponent);
        }

    T operator/ (const Big &other) const
        {
            T rv = this->m_base / other.m_base;
            if (this->m_exponent != other.m_exponent)
                rv *= std::exp(this->m_exponent - other.m_exponent);
            return rv;
        }

    bool is_zero (void) const
        {
            return m_base == T(0);
        }

    bool is_nonzero (void) const
        {
            return m_base != T(0);
        }

    const T & get_base (void) const
        {
            return m_base;
        }

    const exponent_t & get_exponent (void) const
        {
            return m_exponent;
        }

    /**
     * Returns the actual value represented.
     *
     * WARNING: in many cases this is likely to overflow!  Avoiding overflow is
     * the whole point of this class.  But then again, getting the actual value
     * can be useful for testing purposes, and so this method exists.
     */
    T get_value (void) const
        {
            return m_base * std::exp(m_exponent);
        }

private:
    T m_base;
    exponent_t m_exponent;
};

#endif
