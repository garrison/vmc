from libc.string cimport const_char
from libcpp.string cimport string

cdef extern from "vmc-core.hpp":
    cdef cppclass CppHighlevelSimulation "HighlevelSimulation":
        CppHighlevelSimulation(const_char*) except +
        void iterate(int)
        string output()

cdef class HighlevelSimulation(object):
    cdef CppHighlevelSimulation *thisptr

    def __init__(self, input_str):
        cdef unicode input_unicode = unicode(input_str)
        cdef bytes input_bytes = input_unicode.encode('UTF-8')
        cdef char* input_cstr = input_bytes
        self.thisptr = new CppHighlevelSimulation(input_cstr)

    def iterate(self, int sweeps):
        self.thisptr.iterate(sweeps)

    def output(self):
        cdef string std_string = self.thisptr.output()
        # CYTHON-LIMITATION: fixme: casting away constness on the next line
        # seems to be the only way to get it to work, as const_char* doesn't
        # have a "decode" method or subscript operator, which we use on the
        # following line.
        cdef char* c_string = <char*>std_string.c_str()
        return c_string[:std_string.size()].decode('UTF-8', 'strict')

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr
