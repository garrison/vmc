#ifndef _VMC_LW_VECTOR_HPP
#define _VMC_LW_VECTOR_HPP

#ifndef VMC_LW_VECTOR_IS_STD_VECTOR

#include <cstddef>
#include <utility>

#include <boost/array.hpp>

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
    lw_vector (void)
        {
        }

    lw_vector (const lw_vector &lwv)
    : v(lwv.v),
      n(lwv.n)
        {
        }

    lw_vector (unsigned int initial_size, const T &value=T())
    : n(initial_size)
        {
            BOOST_ASSERT(n <= MAX_SIZE);
            v.assign(value);
        }

    lw_vector (const std::initializer_list<T> &l)
        {
            for (const T & x : l)
                this->push_back(x);
        }

    std::size_t size (void) const
        {
            return n;
        }

    void push_back (const T &value)
        {
            ++n;
            BOOST_ASSERT(n <= MAX_SIZE);
            v[n - 1] = value;
        }

    void resize (unsigned int new_size, const T &value=T())
        {
            BOOST_ASSERT(new_size <= MAX_SIZE);
            while (new_size > n) {
                v[n] = value;
                ++n;
            }
            // in case new_size <= n
            n = new_size;
        }

    T & operator[] (std::size_t index)
        {
            BOOST_ASSERT(index < n);
            return v[index];
        }

    const T & operator[] (std::size_t index) const
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
            for (unsigned int i = 0; i < n; ++i) {
                if (v[i] != other.v[i])
                    return v[i] < other.v[i];
            }
            return false;
        }

    typedef const T *const_iterator;

    const_iterator begin (void) const
        {
            return v.data();
        }

    const_iterator end (void) const
        {
            return v.data() + n;
        }

private:
    boost::array<T, MAX_SIZE> v;
    unsigned int n = 0;
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
