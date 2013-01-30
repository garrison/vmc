from pyvmc.includes.libcpp.memory cimport auto_ptr

cdef extern from "RandomNumberGenerator.hpp":
    cdef cppclass CppRandomNumberGenerator "RandomNumberGenerator":
        pass

    bint rng_name_is_valid "RandomNumberGenerator::name_is_valid" (const char*)
    auto_ptr[CppRandomNumberGenerator] create_rng "RandomNumberGenerator::create" (const char*, unsigned long)

cdef class RandomNumberGenerator(object):
    cdef auto_ptr[CppRandomNumberGenerator] autoptr
    cdef unsigned long seed
