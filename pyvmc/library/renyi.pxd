from pyvmc.includes.libcpp.memory cimport unique_ptr, shared_ptr

from pyvmc.core.walk cimport CppWalk
from pyvmc.core.measurement cimport CppBaseMeasurement
from pyvmc.core.estimate cimport CppRealBlockedEstimate, CppComplexBlockedEstimate
from pyvmc.core.wavefunction cimport CppWavefunctionAmplitude
from pyvmc.core.subsystem cimport CppSubsystem

cdef extern from "RenyiModPossibleWalk.hpp":
    cdef cppclass CppRenyiModPossibleWalk "RenyiModPossibleWalk" (CppWalk):
        CppRenyiModPossibleWalk(unique_ptr[CppWavefunctionAmplitude], unique_ptr[CppWavefunctionAmplitude], shared_ptr[CppSubsystem]&)
        CppRenyiModPossibleWalk(unique_ptr[CppWavefunctionAmplitude], shared_ptr[CppSubsystem]&)

cdef extern from "RenyiSignWalk.hpp":
    cdef cppclass CppRenyiSignWalk "RenyiSignWalk" (CppWalk):
        CppRenyiSignWalk(unique_ptr[CppWavefunctionAmplitude], unique_ptr[CppWavefunctionAmplitude], shared_ptr[CppSubsystem]&)
        CppRenyiSignWalk(unique_ptr[CppWavefunctionAmplitude], shared_ptr[CppSubsystem]&)

cdef extern from "RenyiModPossibleMeasurement.hpp":
    cdef cppclass CppRenyiModPossibleMeasurement "RenyiModPossibleMeasurement" (CppBaseMeasurement):
        CppRenyiModPossibleMeasurement()
        CppRealBlockedEstimate& get_estimate()

cdef extern from "RenyiSignMeasurement.hpp":
    cdef cppclass CppRenyiSignMeasurement "RenyiSignMeasurement" (CppBaseMeasurement):
        CppRenyiSignMeasurement()
        CppComplexBlockedEstimate& get_estimate()
