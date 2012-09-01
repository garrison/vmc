from pyvmc.libcpp.memory cimport auto_ptr
from pyvmc.boost.shared_ptr cimport shared_ptr

cdef extern from "Lattice.hpp":
    cdef unsigned int MAX_DIMENSION

    cdef cppclass CppLatticeSite "LatticeSite":
        CppLatticeSite(unsigned int)
        CppLatticeSite(CppLatticeSite&)

        int operator[](int)
        void set_bs_coordinate(int, int)
        int n_dimensions()

        int basis_index

    cdef cppclass lw_vector[T]:
        T& operator[](int)
        int size()
        void push_back(T&)

    ctypedef lw_vector[int] DimensionVector "lw_vector<int, MAX_DIMENSION>"
    ctypedef DimensionVector const_DimensionVector "const lw_vector<int, MAX_DIMENSION>"

    ctypedef lw_vector[int] UDimensionVector "lw_vector<unsigned int, MAX_DIMENSION>"
    ctypedef UDimensionVector const_UDimensionVector "const lw_vector<unsigned int, MAX_DIMENSION>"

    cdef cppclass CppLattice "Lattice":
        CppLattice(DimensionVector, int)

        int total_sites()
        int n_dimensions()
        CppLatticeSite site_from_index(unsigned int)
        unsigned int site_to_index(CppLatticeSite)
        bint site_is_valid(CppLatticeSite)

        DimensionVector dimensions
        int basis_indices

cdef class LatticeSite(object):
    cdef auto_ptr[CppLatticeSite] autoptr

cdef class Lattice(object):
    cdef shared_ptr[CppLattice] sharedptr
