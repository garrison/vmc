from cython.operator cimport dereference as deref
from libcpp.list cimport list as stdlist
from pyvmc.includes.libcpp.memory cimport auto_ptr
from pyvmc.includes.boost.shared_ptr cimport shared_ptr

import logging
import collections

from pyvmc.core.rng cimport RandomNumberGenerator, CppRandomNumberGenerator
from pyvmc.core.lattice cimport Lattice, CppLattice
from pyvmc.core.walk cimport Walk, CppWalk
from pyvmc.core.walk import WalkPlan
from pyvmc.core.measurement cimport BaseMeasurement, CppBaseMeasurement
from pyvmc.utils.resource_logging import log_rusage

cdef extern from "MetropolisSimulation.hpp":
    cdef cppclass CppMetropolisSimulation "MetropolisSimulation":
        CppMetropolisSimulation(auto_ptr[CppWalk], stdlist[shared_ptr[CppBaseMeasurement]], unsigned int, auto_ptr[CppRandomNumberGenerator]&) nogil
        void iterate(unsigned int) nogil
        unsigned int steps_completed()
        unsigned int steps_accepted()
        unsigned int steps_fully_rejected()
        void reset_measurement_estimates()

logger = logging.getLogger(__name__)

cdef class MetropolisSimulation(object):
    cdef auto_ptr[CppMetropolisSimulation] autoptr
    cdef object _walk_plan

    def __init__(self, walk_plan not None, Lattice lattice not None, measurements, unsigned int equilibrium_steps, RandomNumberGenerator rng not None):
        """keep in mind that the rng passed can no longer be used for other things afterwards"""
        assert rng.is_good()

        assert isinstance(walk_plan, WalkPlan)
        self._walk_plan = walk_plan
        cdef Walk walk = walk_plan.create_walk(rng)
        assert walk.autoptr.get() is not NULL
        cdef auto_ptr[CppWalk] walk_autoptr = walk.autoptr

        cdef stdlist[shared_ptr[CppBaseMeasurement]] measurement_list
        cdef BaseMeasurement measurement_
        assert isinstance(measurements, collections.Sequence)
        for measurement in measurements:
            measurement_ = measurement
            if not measurement_.sharedptr.get().is_valid_walk(deref(walk_autoptr)):
                raise ValueError("invalid walk/measurement/wavefunction combination")
            measurement_list.push_back(measurement_.sharedptr)

        with log_rusage(logger, "Equilibrated walk using {} steps.".format(equilibrium_steps)):
            with nogil:
                self.autoptr.reset(new CppMetropolisSimulation(walk_autoptr, measurement_list, equilibrium_steps, rng.autoptr))

        logger.info("Now prepared to consider %d different measurements.", len(measurements))

    def iterate(self, unsigned int sweeps):
        with log_rusage(logger, "Performed {} sweeps on walk.".format(sweeps)):
            with nogil:
                self.autoptr.get().iterate(sweeps)

    def reset_measurement_estimates(self):
        self.autoptr.get().reset_measurement_estimates()

    property walk_plan:
        def __get__(self):
            return self._walk_plan

    property steps_completed:
        def __get__(self):
            return self.autoptr.get().steps_completed()

    property steps_accepted:
        def __get__(self):
            return self.autoptr.get().steps_accepted()

    property steps_fully_rejected:
        def __get__(self):
            return self.autoptr.get().steps_fully_rejected()

    property run_information:
        def __get__(self):
            return RunInformation()

cdef extern from "RunInformation.hpp" namespace "RunInformation":
    const char *compiler
    const char *boost_version
    const char *eigen_version

    int precision_digits "RunInformation::Precision::digits"
    int precision_min_exponent "RunInformation::Precision::min_exponent"
    int precision_max_exponent "RunInformation::Precision::max_exponent"

cdef class RunInformation(object):
    property compiler:
        def __get__(self):
            return compiler

    property boost_version:
        def __get__(self):
            return boost_version

    property eigen_version:
        def __get__(self):
            return eigen_version

    property precision:
        def __get__(self):
            return PrecisionInformation()

cdef class PrecisionInformation(object):
    property digits:
        def __get__(self):
            return precision_digits

    property min_exponent:
        def __get__(self):
            return precision_min_exponent

    property max_exponent:
        def __get__(self):
            return precision_max_exponent
