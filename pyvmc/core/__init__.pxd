cdef extern from "vmc-typedefs.hpp":
    cdef cppclass complex_t:
        float real()
        float imag()
