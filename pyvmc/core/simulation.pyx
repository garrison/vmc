from cython.operator cimport dereference as deref
from libc.string cimport const_char
from libcpp.string cimport string
from libcpp.list cimport list as stdlist

import collections

from pyvmc.boost.shared_ptr cimport shared_ptr
from pyvmc.core.lattice cimport Lattice, CppLattice
from pyvmc.core.measurement cimport BaseMeasurement, CppBaseMeasurement

cdef extern from "MetropolisSimulation.hpp":
    cdef cppclass CppHighlevelSimulation "MetropolisSimulation":
        void iterate(int) nogil
        int steps_completed()
        int steps_accepted()
        int steps_fully_rejected()

cdef extern from "vmc-core.hpp":
        CppHighlevelSimulation* create_simulation(const_char*, shared_ptr[CppLattice], stdlist[shared_ptr[CppBaseMeasurement]], int) except +
        string simulation_output(CppHighlevelSimulation*)

cdef class HighlevelSimulation(object):
    cdef CppHighlevelSimulation *thisptr

    def __init__(self, input_str, Lattice lattice not None, measurements, int equilibrium_steps):
        cdef unicode input_unicode = unicode(input_str)
        cdef bytes input_bytes = input_unicode.encode('UTF-8')
        cdef char* input_cstr = input_bytes
        cdef stdlist[shared_ptr[CppBaseMeasurement]] measurement_list
        cdef BaseMeasurement measurement_
        assert isinstance(measurements, collections.Sequence)
        for measurement in measurements:
            measurement_ = measurement
            #if not measurement_.is_valid_walk(xxx):
            #    raise ValueError("invalid walk/measurement combination")
            measurement_list.push_back(deref(measurement_.sharedptr))
        if self.thisptr is not NULL:
            del self.thisptr
        self.thisptr = create_simulation(input_cstr, deref(lattice.sharedptr), measurement_list, equilibrium_steps)

    def __dealloc__(self):
        del self.thisptr

    def iterate(self, int sweeps):
        self.thisptr.iterate(sweeps)

    def output(self):
        cdef string std_string = simulation_output(self.thisptr)
        # CYTHON-LIMITATION: fixme: casting away constness on the next line
        # seems to be the only way to get it to work, as const_char* doesn't
        # have a "decode" method or subscript operator, which we use on the
        # following line.
        cdef char* c_string = <char*>std_string.c_str()
        return c_string[:std_string.size()].decode('UTF-8', 'strict')

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    property steps_completed:
        def __get__(self):
            return self.thisptr.steps_completed()

    property steps_accepted:
        def __get__(self):
            return self.thisptr.steps_accepted()

    property steps_fully_rejected:
        def __get__(self):
            return self.thisptr.steps_fully_rejected()
