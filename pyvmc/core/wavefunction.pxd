from libc.string cimport const_char
from pyvmc.includes.libcpp.memory cimport auto_ptr
from pyvmc.includes.boost.shared_ptr cimport shared_ptr

from pyvmc.core.rng cimport RandomNumberGenerator, CppRandomNumberGenerator
from pyvmc.core.lattice cimport Lattice, CppLattice

cdef extern from "Wavefunction.hpp":
    cdef cppclass CppWavefunctionAmplitude "Wavefunction::Amplitude":
        shared_ptr[CppWavefunctionAmplitude] clone()

cdef extern from "vmc-core.hpp":
    shared_ptr[CppWavefunctionAmplitude] create_wfa_from_json(const_char*, shared_ptr[CppLattice], auto_ptr[CppRandomNumberGenerator]&) except +

cdef shared_ptr[CppWavefunctionAmplitude] create_wfa(wf)
