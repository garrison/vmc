cdef extern from "vmc-typedefs.hpp":
    cdef cppclass complex_t:
        complex_t()
        complex_t(double, double)
        float real()
        float imag()
