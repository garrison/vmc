import abc
import numbers
import collections

from cython.operator cimport dereference as deref

import numpy

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.subsystem cimport Subsystem
from pyvmc.core.lattice cimport Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.walk import WalkPlan, StandardWalkPlan
from pyvmc.utils.immutable import Immutable

class MeasurementPlan(Immutable):
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

cdef class BaseMeasurement(object):
    pass

class SiteHop(Immutable):
    __slots__ = ("source", "destination", "species")

    def init_validate(self, source, destination, species):
        assert isinstance(source, LatticeSite)
        assert isinstance(destination, LatticeSite)
        assert isinstance(species, numbers.Integral) and species >= 0
        return (source, destination, species)

    def is_valid_for(self, wavefunction):
        assert isinstance(wavefunction, Wavefunction)
        return (self.source in wavefunction.lattice and
                self.destination in wavefunction.lattice and
                self.species < wavefunction.N_species)

class OperatorMeasurementPlan(MeasurementPlan):
    __slots__ = ("walk", "hops", "sum", "boundary_conditions", "steps_per_measurement")

    def __init__(self, wavefunction, hops, sum, boundary_conditions, steps_per_measurement=1000):
        walk = StandardWalkPlan(wavefunction)
        assert isinstance(hops, collections.Sequence)
        hops = tuple(sorted(hops, key=lambda hop: (hop.species, hop.source, hop.destination)))
        assert all([isinstance(hop, SiteHop) and hop.is_valid_for(wavefunction) for hop in hops])
        assert isinstance(sum, bool)
        if not sum and boundary_conditions is not None:
            from warnings import warn
            warn("boundary_conditions are needlessly given, as sum==False", RuntimeWarning)
            boundary_conditions = None
        if boundary_conditions is not None:
            assert valid_boundary_conditions(boundary_conditions, len(wavefunction.lattice.dimensions))
            boundary_conditions = tuple(boundary_conditions)
        super(OperatorMeasurementPlan, self).__init__(walk, hops, sum, boundary_conditions, steps_per_measurement)

    def to_json(self):
        lattice = self.walk.wavefunction.lattice
        hops_list = [{
            "source": lattice.index(hop.source),
            "destination": lattice.index(hop.destination),
            "species": hop.species,
        } for hop in self.hops]
        return {
            "type": "operator",
            "hops": hops_list,
            "sum": self.sum,
            "boundary-conditions": self.boundary_conditions,
            "steps-per-measurement": self.steps_per_measurement,
        }

    def to_measurement(self):
        return OperatorMeasurement(self.steps_per_measurement, self.hops, self.sum, self.boundary_conditions, self.walk.wavefunction.lattice)

cdef class OperatorMeasurement(BaseMeasurement):
    def __init__(self, int steps_per_measurement, hops, sum, bcs, Lattice lattice not None):
        cdef CppBoundaryConditions *cppbcs = NULL
        cdef CppBoundaryConditions cppbcs_

        cdef vector[CppSiteHop] hopv
        cdef LatticeSite src, dest
        for hop in hops:
            assert isinstance(hop, SiteHop)
            src = hop.source
            dest = hop.destination
            hopv.push_back(CppSiteHop(src.cpp, dest.cpp, hop.species))
        cdef CppParticleOperator *operator
        operator = new CppParticleOperator(hopv, lattice.sharedptr)
        try:
            if bcs:
                for bc in bcs:
                    # NOTE: we store the fraction's inverse in python vs c++ code
                    cppbcs_.push_back(CppBoundaryCondition(boost_rational[int](bc.denominator, bc.numerator)))
                cppbcs = &cppbcs_
            self.sharedptr.reset(new CppOperatorMeasurement(steps_per_measurement, deref(operator), sum, cppbcs))
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
        return {
            "type": "subsystem-occupation-number-probability",
            "subsystem": self.subsystem.to_json(),
            "steps-per-measurement": self.steps_per_measurement,
        }

    def to_measurement(self):
        return SubsystemOccupationNumberProbabilityMeasurement(self.steps_per_measurement, self.subsystem)

cdef class SubsystemOccupationNumberProbabilityMeasurement(BaseMeasurement):
    def __init__(self, int steps_per_measurement, Subsystem subsystem not None):
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
