from cython.operator cimport dereference as deref
from libcpp.list cimport list as stdlist
from pyvmc.includes.libcpp.memory cimport auto_ptr
from pyvmc.includes.boost.shared_ptr cimport shared_ptr

import collections

from pyvmc.core.rng cimport RandomNumberGenerator, CppRandomNumberGenerator
from pyvmc.core.lattice cimport Lattice, CppLattice
from pyvmc.core.walk cimport Walk, CppWalk
from pyvmc.core.measurement cimport BaseMeasurement, CppBaseMeasurement

cdef extern from "MetropolisSimulation.hpp":
    cdef cppclass CppMetropolisSimulation "MetropolisSimulation":
        CppMetropolisSimulation(auto_ptr[CppWalk], stdlist[shared_ptr[CppBaseMeasurement]], unsigned int, auto_ptr[CppRandomNumberGenerator]&) nogil
        void iterate(unsigned int) nogil
        unsigned int steps_completed()
        unsigned int steps_accepted()
        unsigned int steps_fully_rejected()
        void reset_measurement_estimates()

cdef class MetropolisSimulation(object):
    cdef auto_ptr[CppMetropolisSimulation] autoptr

    def __init__(self, Walk walk not None, Lattice lattice not None, measurements, unsigned int equilibrium_steps):
        if walk.autoptr.get() is NULL:
            raise RuntimeError("Walk's auto_ptr is null.  It cannot be recycled.")
        cdef auto_ptr[CppWalk] walk_autoptr = walk.autoptr

        cdef stdlist[shared_ptr[CppBaseMeasurement]] measurement_list
        cdef BaseMeasurement measurement_
        assert isinstance(measurements, collections.Sequence)
        for measurement in measurements:
            measurement_ = measurement
            if not measurement_.sharedptr.get().is_valid_walk(deref(walk_autoptr)):
                raise ValueError("invalid walk/measurement/wavefunction combination")
            measurement_list.push_back(measurement_.sharedptr)

        rng = RandomNumberGenerator()
        cdef auto_ptr[CppRandomNumberGenerator] rng_autoptr = rng.autoptr

        with nogil:
            self.autoptr.reset(new CppMetropolisSimulation(walk_autoptr, measurement_list, equilibrium_steps, rng_autoptr))

    def iterate(self, unsigned int sweeps):
        with nogil:
            self.autoptr.get().iterate(sweeps)

    def reset_measurement_estimates(self):
        self.autoptr.get().reset_measurement_estimates()

    property steps_completed:
        def __get__(self):
            return self.autoptr.get().steps_completed()

    property steps_accepted:
        def __get__(self):
            return self.autoptr.get().steps_accepted()

    property steps_fully_rejected:
        def __get__(self):
            return self.autoptr.get().steps_fully_rejected()
