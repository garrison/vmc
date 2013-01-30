from pyvmc.includes.boost.shared_ptr cimport shared_ptr

from pyvmc.core.lattice cimport CppLattice, UDimensionVector

cdef extern from "Subsystem.hpp":
    cdef unsigned int MAX_DIMENSION

    cdef cppclass CppSubsystem "Subsystem":
        bint position_is_within(unsigned int site_index, CppLattice &lattice)
        bint lattice_makes_sense(CppLattice &lattice)

cdef extern from "SimpleSubsystem.hpp":
    cdef cppclass CppSimpleSubsystem "SimpleSubsystem" (CppSubsystem):
        CppSimpleSubsystem(UDimensionVector)

        UDimensionVector subsystem_length

cdef class Subsystem(object):
    cdef object lattice_
    cdef shared_ptr[CppSubsystem] sharedptr
