cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T*)
        void reset()
        void reset(T*)
        T& operator*()
        T* get()

cdef extern from "Lattice.hpp":
    cdef unsigned int MAX_DIMENSION

    cdef cppclass CppLatticeSite "LatticeSite":
        CppLatticeSite(unsigned int)
        CppLatticeSite(CppLatticeSite&)

        int operator[](int)
        void set_bs_coordinate(int, int)
        int n_dimensions()

        int basis_index

    cppclass DimensionVector:
        int operator[](int)
        void push_back(int)

    cdef cppclass CppLattice "Lattice":
        CppLattice(DimensionVector, int)

        int total_sites()
        int n_dimensions()
        CppLatticeSite site_from_index(unsigned int)
        unsigned int site_to_index(CppLatticeSite)
        int site_is_valid(CppLatticeSite)

        DimensionVector dimensions
        int basis_indices

cdef class LatticeSite(object):
    cdef CppLatticeSite *thisptr

cdef class Lattice(object):
    cdef shared_ptr[CppLattice] *sharedptr
