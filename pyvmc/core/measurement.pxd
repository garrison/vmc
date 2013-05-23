from pyvmc.includes.boost.shared_ptr cimport shared_ptr

from pyvmc.core.walk cimport CppWalk

cdef extern from "Measurement.hpp":
    cdef cppclass CppBaseMeasurement "BaseMeasurement<probability_t>":
        bint is_valid_walk(CppWalk&)

cdef class BaseMeasurement(object):
    cdef shared_ptr[CppBaseMeasurement] sharedptr
