#ifndef _VMC_LW_VECTOR_HPP
#define _VMC_LW_VECTOR_HPP

#ifndef VMC_LW_VECTOR_IS_STD_VECTOR

#include <cstddef>
#include <array>
#include <utility>

#include <boost/assert.hpp>

/**
 * Lightweight vector (no dynamic allocation) with a maximum size
 *
 * This is meant to act like a standard STL container, but we only implement
 * things as they are needed.
 */
template <typename T, unsigned int MAX_SIZE>
class lw_vector
{
public:
    typedef T value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef pointer iterator;
    typedef const_pointer const_iterator;
    typedef unsigned int size_type;
    typedef std::ptrdiff_t difference_type;

    lw_vector (void)
        {
        }

    lw_vector (const lw_vector &lwv)
    : v(lwv.v),
      n(lwv.n)
        {
        }

    lw_vector (size_type initial_size, const T &value=T())
    : n(initial_size)
        {
            BOOST_ASSERT(n <= MAX_SIZE);
            v.fill(value);
        }

    lw_vector (const std::initializer_list<T> &l)
        {
            for (const T & x : l)
                this->push_back(x);
        }

    size_type size (void) const
        {
            return n;
        }

    size_type max_size (void) const
        {
            return MAX_SIZE;
        }

    bool empty (void) const
        {
            return n == 0;
        }

    void push_back (const T &value)
        {
            ++n;
            BOOST_ASSERT(n <= MAX_SIZE);
            v[n - 1] = value;
        }

    void resize (size_type new_size, const T &value=T())
        {
            BOOST_ASSERT(new_size <= MAX_SIZE);
            while (new_size > n) {
                v[n] = value;
                ++n;
            }
            // in case new_size <= n
            n = new_size;
        }

    reference operator[] (size_type index)
        {
            BOOST_ASSERT(index < n);
            return v[index];
        }

    const_reference operator[] (size_type index) const
        {
            BOOST_ASSERT(index < n);
            return v[index];
        }

    bool operator== (const lw_vector<T, MAX_SIZE> &other) const
        {
            return (n == other.n
                    && v == other.v);
        }

    bool operator!= (const lw_vector<T, MAX_SIZE> &other) const
        {
            return (n != other.n
                    || v != other.v);
        }

    bool operator< (const lw_vector<T, MAX_SIZE> &other) const
        {
            if (n != other.n)
                return n < other.n;
            for (size_type i = 0; i < n; ++i) {
                if (v[i] != other.v[i])
                    return v[i] < other.v[i];
            }
            return false;
        }

    iterator begin (void)
        {
            return v.data();
        }

    const_iterator begin (void) const
        {
            return v.data();
        }

    iterator end (void)
        {
            return v.data() + n;
        }

    const_iterator end (void) const
        {
            return v.data() + n;
        }

    const_iterator cbegin (void) const
        {
            return begin();
        }

    const_iterator cend (void) const
        {
            return end();
        }

private:
    std::array<T, MAX_SIZE> v;
    size_type n = 0;
};

#else // VMC_LW_VECTOR_IS_STD_VECTOR

#include <vector>

template <typename T, unsigned int MAX_SIZE>
class lw_vector : public std::vector<T>
{
public:
    lw_vector ()
        : std::vector<T>()
        {
        }

    lw_vector (std::size_t size)
        : std::vector<T>(size)
        {
        }
};

#endif // VMC_LW_VECTOR_IS_STD_VECTOR

#endif
