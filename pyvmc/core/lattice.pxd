from pyvmc.core cimport lw_vector
from pyvmc.includes.libcpp.memory cimport unique_ptr, shared_ptr

cdef extern from "Lattice.hpp":
    cdef unsigned int MAX_DIMENSION

    cdef cppclass CppLatticeSite "LatticeSite":
        CppLatticeSite()
        CppLatticeSite(unsigned int)
        CppLatticeSite(const CppLatticeSite&)

        int operator[](int)
        void set_n_dimensions(unsigned int)
        void set_bs_coordinate(size_t, int)
        unsigned int n_dimensions()

        int basis_index

        bint operator==(const CppLatticeSite&)
        bint operator!=(const CppLatticeSite&)
        bint operator<(const CppLatticeSite&)

    ctypedef lw_vector[int] CppLatticeDimensionVector "Lattice::DimensionVector"

    cdef cppclass CppLattice "Lattice":
        CppLattice(const CppLatticeDimensionVector&, int)

        unsigned int total_sites()
        unsigned int total_bravais_sites()
        unsigned int n_dimensions()
        CppLatticeSite operator[](unsigned int)
        unsigned int index(const CppLatticeSite&)
        bint site_is_valid(const CppLatticeSite&)

        const CppLatticeDimensionVector dimensions
        const int basis_indices

cdef class LatticeSite(object):
    cdef CppLatticeSite cpp

cdef class Lattice(object):
    cdef shared_ptr[CppLattice] sharedptr
