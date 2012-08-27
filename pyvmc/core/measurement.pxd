from pyvmc.boost.shared_ptr cimport shared_ptr
from pyvmc.core.subsystem cimport CppSubsystem

cdef extern from "Measurement.hpp":
    cdef cppclass CppBaseMeasurement "BaseMeasurement":
        pass

cdef class BaseMeasurement(object):
    cdef shared_ptr[CppBaseMeasurement] *sharedptr

cdef extern from "SubsystemOccupationNumberProbabilityMeasurement.hpp":
    cdef cppclass CppSubsystemOccupationNumberProbabilityMeasurement "SubsystemOccupationNumberProbabilityMeasurement" (CppBaseMeasurement):
        CppSubsystemOccupationNumberProbabilityMeasurement(int, shared_ptr[CppSubsystem])
