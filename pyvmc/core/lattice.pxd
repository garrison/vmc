from pyvmc.includes.libcpp.memory cimport auto_ptr
from pyvmc.includes.boost.shared_ptr cimport shared_ptr

cdef extern from "Lattice.hpp":
    cdef unsigned int MAX_DIMENSION

    cdef cppclass CppLatticeSite "LatticeSite":
        CppLatticeSite()
        CppLatticeSite(unsigned int)
        CppLatticeSite(CppLatticeSite&)

        int operator[](int)
        void set_n_dimensions(unsigned int)
        void set_bs_coordinate(size_t, int)
        unsigned int n_dimensions()

        int basis_index

    cdef cppclass lw_vector[T]:
        T& operator[](size_t)
        size_t size()
        void push_back(T&)

    ctypedef lw_vector[int] DimensionVector "lw_vector<int, MAX_DIMENSION>"
    ctypedef DimensionVector const_DimensionVector "const lw_vector<int, MAX_DIMENSION>"

    ctypedef lw_vector[int] UDimensionVector "lw_vector<unsigned int, MAX_DIMENSION>"
    ctypedef UDimensionVector const_UDimensionVector "const lw_vector<unsigned int, MAX_DIMENSION>"

    cdef cppclass CppLattice "Lattice":
        CppLattice(DimensionVector, int)

        unsigned int total_sites()
        unsigned int n_dimensions()
        CppLatticeSite site_from_index(unsigned int)
        unsigned int site_to_index(CppLatticeSite)
        bint site_is_valid(CppLatticeSite)

        DimensionVector dimensions
        int basis_indices

cdef class LatticeSite(object):
    cdef CppLatticeSite cpp

cdef class Lattice(object):
    cdef shared_ptr[CppLattice] sharedptr
