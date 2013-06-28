from libcpp.vector cimport vector
from pyvmc.includes.libcpp.memory cimport shared_ptr

from pyvmc.core.rng cimport RandomNumberGenerator, CppRandomNumberGenerator
from pyvmc.core.lattice cimport Lattice, CppLattice
from pyvmc.core.orbitals cimport CppOrbitalDefinitions, const_CppOrbitalDefinitions

cdef extern from "Wavefunction.hpp":
    cdef cppclass CppWavefunctionAmplitude "Wavefunction<amplitude_t>::Amplitude":
        shared_ptr[CppWavefunctionAmplitude] clone()

    cdef cppclass CppWavefunction "Wavefunction<amplitude_t>":
        shared_ptr[CppWavefunctionAmplitude] create_nonzero_wavefunctionamplitude(shared_ptr[CppWavefunction]&, CppRandomNumberGenerator&) except +
        shared_ptr[CppWavefunctionAmplitude] create_nonzero_wavefunctionamplitude(shared_ptr[CppWavefunction]&, CppRandomNumberGenerator&, unsigned int) except +

cdef shared_ptr[CppWavefunctionAmplitude] create_nonzero_wfa(wf, RandomNumberGenerator rng) except *

cdef extern from "FreeFermionWavefunction.hpp":
    cdef cppclass CppFreeFermionWavefunction "FreeFermionWavefunction<amplitude_t>" (CppWavefunction):
        CppFreeFermionWavefunction(vector[shared_ptr[const_CppOrbitalDefinitions]]&, shared_ptr[CppJastrowFactor]&)

cdef class WavefunctionWrapper(object):
    cdef shared_ptr[CppWavefunction] sharedptr

cdef extern from "JastrowFactor.hpp":
    cdef cppclass CppJastrowFactor "JastrowFactor<amplitude_t>":
        pass

cdef extern from "NoDoubleOccupancyProjector.hpp":
    cdef cppclass CppNoDoubleOccupancyProjector "NoDoubleOccupancyProjector<amplitude_t>" (CppJastrowFactor):
        CppNoDoubleOccupancyProjector()

cdef extern from "JordanWignerJastrowFactor.hpp":
    cdef cppclass CppJordanWignerJastrowFactor "JordanWignerJastrowFactor<amplitude_t>" (CppJastrowFactor):
        CppJordanWignerJastrowFactor()

cdef extern from "TwoBodyJastrowFactor.hpp":
    cdef cppclass CppTwoBodyJastrowMatrix "Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic>":
        CppTwoBodyJastrowMatrix(unsigned int, unsigned int)
        CppTwoBodyJastrowMatrix(const CppTwoBodyJastrowMatrix&)

    cdef cppclass CppTwoBodyJastrowFactor "TwoBodyJastrowFactor<amplitude_t>" (CppJastrowFactor):
        CppTwoBodyJastrowFactor(const CppTwoBodyJastrowMatrix&)

    cdef void set_matrix_coeff "TwoBodyJastrowFactor<amplitude_t>::set_matrix_coeff" (CppTwoBodyJastrowMatrix&, unsigned int, unsigned int, double)

cdef class JastrowFactor(object):
    cdef shared_ptr[CppJastrowFactor] sharedptr
