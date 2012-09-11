from libc.string cimport const_char
from pyvmc.includes.libcpp.memory cimport auto_ptr

cdef extern from "RandomNumberGenerator.hpp":
    cdef cppclass CppRandomNumberGenerator "RandomNumberGenerator":
        pass

    bint rng_name_is_valid "RandomNumberGenerator::name_is_valid" (const_char*)
    auto_ptr[CppRandomNumberGenerator] create_rng "RandomNumberGenerator::create" (const_char*, unsigned long)

cdef class RandomNumberGenerator(object):
    cdef auto_ptr[CppRandomNumberGenerator] autoptr
    cdef unsigned long seed
