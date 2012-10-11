cdef extern from "<boost/rational.hpp>" namespace "boost":
    cdef cppclass rational[T]:
        rational() nogil
        rational(T) nogil
        rational(T, T) nogil
        void assign(T, T) nogil
        T numerator() nogil
        T denominator() nogil
