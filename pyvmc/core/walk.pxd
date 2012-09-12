from pyvmc.includes.libcpp.memory cimport auto_ptr
from pyvmc.includes.boost.shared_ptr cimport shared_ptr

from pyvmc.core.wavefunction cimport CppWavefunctionAmplitude

cdef extern from "Walk.hpp":
    cdef cppclass CppWalk "Walk":
        pass

cdef extern from "StandardWalk.hpp":
    cdef cppclass CppStandardWalk "StandardWalk" (CppWalk):
        CppStandardWalk(shared_ptr[CppWavefunctionAmplitude]&)

cdef class Walk(object):
    # since it's stored as an autoptr, be careful to check that it's not null
    cdef auto_ptr[CppWalk] autoptr
