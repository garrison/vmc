from pyvmc.includes.libcpp.memory cimport unique_ptr

cdef extern from "RandomNumberGenerator.hpp":
    cdef cppclass CppRandomNumberGenerator "RandomNumberGenerator":
        pass

    bint rng_name_is_valid "RandomNumberGenerator::name_is_valid" (const char*)
    unique_ptr[CppRandomNumberGenerator] create_rng "RandomNumberGenerator::create" (const char*, unsigned long)

cdef class RandomNumberGenerator(object):
    cdef unique_ptr[CppRandomNumberGenerator] autoptr
    cdef str _name
    cdef unsigned long _seed
