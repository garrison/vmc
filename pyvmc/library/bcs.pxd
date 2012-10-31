from libcpp.vector cimport vector
from pyvmc.includes.boost.shared_ptr cimport shared_ptr
from pyvmc.core cimport complex_t

from pyvmc.core.wavefunction cimport CppWavefunction, WavefunctionWrapper
from pyvmc.core.lattice cimport CppLattice, Lattice

cdef extern from "BCSWavefunction.hpp":
    cdef cppclass CppBCSWavefunction "BCSWavefunction" (CppWavefunction):
        CppBCSWavefunction(shared_ptr[CppLattice]&, vector[complex_t]&, unsigned int)
