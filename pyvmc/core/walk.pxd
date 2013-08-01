from pyvmc.includes.libcpp.memory cimport unique_ptr

from pyvmc.core.wavefunction cimport CppWavefunctionAmplitude

cdef extern from "Walk.hpp":
    cdef cppclass CppWalk "Walk<probability_t>":
        pass

cdef extern from "StandardWalk.hpp":
    cdef cppclass CppStandardWalk "StandardWalk<amplitude_t>" (CppWalk):
        CppStandardWalk(unique_ptr[CppWavefunctionAmplitude])

cdef class Walk(object):
    # since it's stored as a unique_ptr, be careful to check that it's not null
    cdef unique_ptr[CppWalk] uniqueptr
