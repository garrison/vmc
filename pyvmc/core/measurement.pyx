import abc
import collections

from cython.operator cimport dereference as deref
from pyvmc.includes.libcpp.memory cimport auto_ptr

import numpy

from pyvmc.core cimport complex_t
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.boundary_conditions cimport BoundaryCondition
from pyvmc.core.subsystem cimport Subsystem
from pyvmc.core.lattice cimport Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.walk import WalkPlan, StandardWalkPlan
from pyvmc.core.operator import SiteHop, BasicOperator, Operator
from pyvmc.core.estimate cimport Estimate_from_CppIntegerBinnedEstimate, Estimate_from_CppComplexBinnedEstimate
from pyvmc.utils.immutable import Immutable, ImmutableMetaclass

class MeasurementPlan(Immutable):
    """base class, for both actual measurements and composite measurements"""

    __slots__ = ()

    @abc.abstractmethod
    def get_measurement_plans(self):
        "should return a set of BasicMeasurementPlan's"
        raise NotImplementedError

class CompositeMeasurementPlan(MeasurementPlan):
    __slots__ = ()

__basic_measurement_plan_registry = {}

class BasicMeasurementPlanMetaclass(ImmutableMetaclass):
    def __init__(cls, name, bases, dct):
        assert name not in __basic_measurement_plan_registry
        __basic_measurement_plan_registry[name] = cls
        super(BasicMeasurementPlanMetaclass, cls).__init__(name, bases, dct)

class BasicMeasurementPlan(MeasurementPlan):
    """base class for fundamental measurements implemented in VMC"""

    __slots__ = ("walk",)

    __metaclass__ = BasicMeasurementPlanMetaclass

    def init_validate(self, walk, *args, **kwargs):
        assert isinstance(walk, WalkPlan)
        return super(BasicMeasurementPlan, self).init_validate(walk, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

    #@abc.abstractmethod
    def from_json(self, json_object):
        raise NotImplementedError

    @abc.abstractmethod
    def to_measurement(self):
        raise NotImplementedError

    def get_measurement_plans(self):
        return {self}

#    def get_result(self, universe):
#        return universe[self].get_result()

# fixme: we could move everything below to pyvmc.measurements

class BasicOperatorMeasurementPlan(BasicMeasurementPlan):
    __slots__ = ("walk", "operator", "steps_per_measurement")

    def __init__(self, wavefunction, operator, steps_per_measurement=1000):
        walk = StandardWalkPlan(wavefunction)
        assert isinstance(operator, BasicOperator)
        assert all([hop.is_valid_for(wavefunction) for hop in operator.hops])
        if operator.boundary_conditions is not None:
            assert valid_boundary_conditions(operator.boundary_conditions, len(wavefunction.lattice.dimensions))
        super(BasicOperatorMeasurementPlan, self).__init__(walk, operator, steps_per_measurement)

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
            ("operator", self.operator.to_json()),
            ("steps-per-measurement", self.steps_per_measurement),
        ])

    def to_measurement(self):
        return BasicOperatorMeasurement(self.steps_per_measurement, self.operator, self.walk.wavefunction.lattice)

cdef class BasicOperatorMeasurement(BaseMeasurement):
    """A VMC measurement which works on BasicOperator's"""

    def __init__(self, unsigned int steps_per_measurement, operator_, Lattice lattice not None):
        assert isinstance(operator_, BasicOperator)

        cdef CppBoundaryConditions cppbcs

        cdef vector[CppSiteHop] hopv
        cdef LatticeSite src, dest
        for hop in operator_.hops:
            assert isinstance(hop, SiteHop)
            src = hop.source
            dest = hop.destination
            hopv.push_back(CppSiteHop(src.cpp, dest.cpp, hop.species))
        cdef auto_ptr[CppParticleOperator] operator
        operator.reset(new CppParticleOperator(hopv, lattice.sharedptr))
        if operator_.boundary_conditions is not None:
            for bc in operator_.boundary_conditions:
                cppbcs.push_back((<BoundaryCondition>bc).cpp)
        self.sharedptr.reset(new CppOperatorMeasurement(steps_per_measurement, deref(operator), cppbcs))

    def get_estimate(self, key=None):
        if key is not None:
            raise KeyError
        return Estimate_from_CppComplexBinnedEstimate((<CppOperatorMeasurement*>self.sharedptr.get()).get_estimate())

    def get_estimates(self):
        return {None: self.get_estimate()}

class SubsystemOccupationProbabilityMeasurementPlan(BasicMeasurementPlan):
    __slots__ = ("walk", "subsystem", "steps_per_measurement")

    def __init__(self, wavefunction, subsystem, steps_per_measurement=100):
        walk = StandardWalkPlan(wavefunction)
        super(SubsystemOccupationProbabilityMeasurementPlan, self).__init__(walk, subsystem, steps_per_measurement)

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
            ("subsystem", self.subsystem.to_json()),
            ("steps-per-measurement", self.steps_per_measurement),
        ])

    def to_measurement(self):
        return SubsystemOccupationNumberProbabilityMeasurement(self.steps_per_measurement, self.subsystem)

cdef class SubsystemOccupationNumberProbabilityMeasurement(BaseMeasurement):
    def __init__(self, unsigned int steps_per_measurement, Subsystem subsystem not None):
        self.sharedptr.reset(new CppSubsystemOccupationNumberProbabilityMeasurement(steps_per_measurement, subsystem.sharedptr))

    def get_estimates(self):
        cdef unsigned int i
        # fixme: in cython 0.17 we will be able to do this iteration directly
        bounds = []
        cdef CppOccupationBounds *cppbounds = &(<CppSubsystemOccupationNumberProbabilityMeasurement*>self.sharedptr.get()).get_bounds()
        for i in xrange(cppbounds.size()):
            # the +1 is because we want to range from 0 to N inclusive
            bounds.append(deref(cppbounds)[i] + 1)
        cdef vector[unsigned int] occ
        occ.resize(len(bounds))
        rv = {}
        for occupation in numpy.ndindex(*bounds):
            for i in xrange(len(bounds)):
                occ[i] = occupation[i]
            rv[tuple(occupation)] = Estimate_from_CppIntegerBinnedEstimate((<CppSubsystemOccupationNumberProbabilityMeasurement*>self.sharedptr.get()).get_estimate(occ))
        return rv

    def get_estimate(self, key):
        # fixme: not efficient
        return self.get_estimates()[key]
