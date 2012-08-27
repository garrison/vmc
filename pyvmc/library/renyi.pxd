from pyvmc.core.measurement cimport CppBaseMeasurement, BaseMeasurement

cdef extern from "RenyiModPossibleMeasurement.hpp":
    cdef cppclass CppRenyiModPossibleMeasurement "RenyiModPossibleMeasurement" (CppBaseMeasurement):
        CppRenyiModPossibleMeasurement()

cdef extern from "RenyiSignMeasurement.hpp":
    cdef cppclass CppRenyiSignMeasurement "RenyiSignMeasurement" (CppBaseMeasurement):
        CppRenyiSignMeasurement()
