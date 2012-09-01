cdef extern from "<boost/rational.hpp>" namespace "boost":
    cdef cppclass rational[T]:
        rational()
        rational(T)
        rational(T, T)
        void assign(T, T)
        T numerator()
        T denominator()
