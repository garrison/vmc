cdef extern from "<boost/dynamic_bitset.hpp>":
    # CYTHON-LIMITATION: no support for parameterless class templates
    cdef cppclass dynamic_bitset "boost::dynamic_bitset<>":
        dynamic_bitset() nogil except +
        dynamic_bitset(const dynamic_bitset &) nogil except +
        dynamic_bitset(size_t) nogil except +
        bint& operator[](size_t) nogil
        size_t size() nogil
        void resize(size_t) nogil
