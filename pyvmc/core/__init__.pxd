cdef extern from "vmc-typedefs.hpp":
    cdef cppclass complex_t:
        complex_t()
        complex_t(double, double)
        double real()
        double imag()

cdef extern from "lw_vector.hpp":
    cdef cppclass lw_vector[T]:
        const T& operator[](size_t)
        size_t size()
        void push_back(T&)
