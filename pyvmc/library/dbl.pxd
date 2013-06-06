from pyvmc.includes.libcpp.memory cimport shared_ptr

from pyvmc.core.orbitals cimport CppOrbitalDefinitions
from pyvmc.core.wavefunction cimport CppWavefunction, WavefunctionWrapper

cdef extern from "DBLWavefunction.hpp":
    cdef cppclass CppDBLWavefunction "DBLWavefunction<amplitude_t>" (CppWavefunction):
        CppDBLWavefunction(shared_ptr[CppOrbitalDefinitions]&, shared_ptr[CppOrbitalDefinitions]&)
        CppDBLWavefunction(shared_ptr[CppOrbitalDefinitions]&, shared_ptr[CppOrbitalDefinitions]&, double, double)
