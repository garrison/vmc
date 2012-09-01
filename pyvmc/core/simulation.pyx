from cython.operator cimport dereference as deref
from libc.string cimport const_char
from libcpp.list cimport list as stdlist
from pyvmc.libcpp.memory cimport auto_ptr
from pyvmc.boost.shared_ptr cimport shared_ptr

import collections

from pyvmc.core.lattice cimport Lattice, CppLattice
from pyvmc.core.measurement cimport BaseMeasurement, CppBaseMeasurement

cdef extern from "MetropolisSimulation.hpp":
    cdef cppclass CppMetropolisSimulation "MetropolisSimulation":
        void iterate(int) nogil
        int steps_completed()
        int steps_accepted()
        int steps_fully_rejected()

cdef extern from "vmc-core.hpp":
        CppMetropolisSimulation* create_simulation(const_char*, shared_ptr[CppLattice], stdlist[shared_ptr[CppBaseMeasurement]], int) except +

cdef class MetropolisSimulation(object):
    cdef auto_ptr[CppMetropolisSimulation] autoptr

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
            measurement_list.push_back(measurement_.sharedptr)
        self.autoptr.reset(create_simulation(input_cstr, lattice.sharedptr, measurement_list, equilibrium_steps))

    def iterate(self, int sweeps):
        self.autoptr.get().iterate(sweeps)

    property steps_completed:
        def __get__(self):
            return self.autoptr.get().steps_completed()

    property steps_accepted:
        def __get__(self):
            return self.autoptr.get().steps_accepted()

    property steps_fully_rejected:
        def __get__(self):
            return self.autoptr.get().steps_fully_rejected()
