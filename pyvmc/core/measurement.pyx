import abc
import collections

from cython.operator cimport dereference as deref

import numpy

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.subsystem cimport Subsystem
from pyvmc.core.lattice cimport Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.walk import WalkPlan, StandardWalkPlan
from pyvmc.core.operator import SiteHop, BasicOperator, Operator
from pyvmc.utils.immutable import Immutable

class BaseMeasurementPlan(Immutable):
    """base class, for both actual measurements and composite measurements"""

    __slots__ = ()

    @abc.abstractmethod
    def get_measurement_plans(self):
        "should return a set"
        raise NotImplementedError

class MeasurementPlan(BaseMeasurementPlan):
    """base class for fundamental measurements implemented in VMC"""

    __slots__ = ("walk",)

    def init_validate(self, walk, *args, **kwargs):
        assert isinstance(walk, WalkPlan)
        return super(MeasurementPlan, self).init_validate(walk, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

    @abc.abstractmethod
    def to_measurement(self):
        raise NotImplementedError

    def get_measurement_plans(self):
        return {self}

    def get_result(self, universe):
        return universe[self].get_result()

cdef class BaseMeasurement(object):
    pass

# fixme: we could move everything below to pyvmc.measurements

class BasicOperatorMeasurementPlan(MeasurementPlan):
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
        return OperatorMeasurement(self.steps_per_measurement, self.operator, self.walk.wavefunction.lattice)

cdef class OperatorMeasurement(BaseMeasurement):
    """A VMC measurement which works on BasicOperator's"""

    def __init__(self, unsigned int steps_per_measurement, operator_, Lattice lattice not None):
        assert isinstance(operator_, BasicOperator)

        cdef CppBoundaryConditions *cppbcs = NULL
        cdef CppBoundaryConditions cppbcs_

        cdef vector[CppSiteHop] hopv
        cdef LatticeSite src, dest
        for hop in operator_.hops:
            assert isinstance(hop, SiteHop)
            src = hop.source
            dest = hop.destination
            hopv.push_back(CppSiteHop(src.cpp, dest.cpp, hop.species))
        cdef CppParticleOperator *operator
        operator = new CppParticleOperator(hopv, lattice.sharedptr)
        try:
            if operator_.boundary_conditions:
                for bc in operator_.boundary_conditions:
                    # NOTE: we store the fraction's inverse in python vs c++ code
                    cppbcs_.push_back(CppBoundaryCondition(boost_rational[int](bc.numerator, bc.denominator)))
                cppbcs = &cppbcs_
            self.sharedptr.reset(new CppOperatorMeasurement(steps_per_measurement, deref(operator), operator_.sum, cppbcs))
        finally:
            del operator

    def get_result(self):
        cdef complex_t c = (<CppOperatorMeasurement*>self.sharedptr.get()).get_estimate().get_result()
        return complex(c.real(), c.imag())

class SubsystemOccupationProbabilityMeasurementPlan(MeasurementPlan):
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

    def get_result(self):
        cdef unsigned int i
        # fixme: in cython 0.17 we will be able to do this iteration directly
        bounds = []
        cdef CppOccupationBounds *cppbounds = &(<CppSubsystemOccupationNumberProbabilityMeasurement*>self.sharedptr.get()).get_bounds()
        for i in xrange(cppbounds.size()):
            # the +1 is because we want to range from 0 to N inclusive
            bounds.append(deref(cppbounds)[i] + 1)
        cdef vector[unsigned int] occ
        occ.resize(len(bounds))
        rv = []
        for occupation in numpy.ndindex(*bounds):
            for i in xrange(len(bounds)):
                occ[i] = occupation[i]
            rv.append((tuple(occupation), (<CppSubsystemOccupationNumberProbabilityMeasurement*>self.sharedptr.get()).get_estimate(occ).get_result()))
        return rv
