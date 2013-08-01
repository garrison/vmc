from pyvmc.includes.libcpp.memory cimport shared_ptr
from pyvmc.includes.boost.dynamic_bitset cimport dynamic_bitset

from pyvmc.core cimport lw_vector
from pyvmc.core.lattice cimport CppLattice

cdef extern from "Subsystem.hpp":
    cdef cppclass CppSubsystem "Subsystem":
        bint position_is_within(unsigned int site_index, CppLattice &lattice)
        bint lattice_makes_sense(CppLattice &lattice)

cdef extern from "SimpleSubsystem.hpp":
    cdef unsigned int MAX_DIMENSION

    ctypedef lw_vector[unsigned int] CppSimpleSubsystemDimensionVector "SimpleSubsystem::DimensionVector"

    cdef cppclass CppSimpleSubsystem "SimpleSubsystem" (CppSubsystem):
        CppSimpleSubsystem(CppSimpleSubsystemDimensionVector)

        CppSimpleSubsystemDimensionVector subsystem_length

cdef extern from "CustomSubsystem.hpp":
    cdef cppclass CppCustomSubsystem "CustomSubsystem" (CppSubsystem):
        CppCustomSubsystem(const dynamic_bitset &)

cdef class Subsystem(object):
    cdef object lattice_
    cdef shared_ptr[CppSubsystem] sharedptr
