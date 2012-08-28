from pyvmc.core.measurement cimport CppBaseMeasurement, CppRealBinnedEstimate, CppComplexBinnedEstimate

cdef extern from "RenyiModPossibleMeasurement.hpp":
    cdef cppclass CppRenyiModPossibleMeasurement "RenyiModPossibleMeasurement" (CppBaseMeasurement):
        CppRenyiModPossibleMeasurement()
        CppRealBinnedEstimate& get_estimate()

cdef extern from "RenyiSignMeasurement.hpp":
    cdef cppclass CppRenyiSignMeasurement "RenyiSignMeasurement" (CppBaseMeasurement):
        CppRenyiSignMeasurement()
        CppComplexBinnedEstimate& get_estimate()
