### declarations (everything here could be moved to a .pxd file if we needed)

from pyvmc.includes.boost.shared_ptr cimport shared_ptr

from libcpp.vector cimport vector

from pyvmc.core.lattice cimport CppLattice, CppLatticeSite
from pyvmc.core.boundary_conditions cimport CppBoundaryCondition, CppBoundaryConditions
from pyvmc.core.subsystem cimport CppSubsystem
from pyvmc.core.estimate cimport CppIntegerBinnedEstimate, CppComplexBinnedEstimate
from pyvmc.core.measurement cimport CppBaseMeasurement

cdef extern from "SubsystemOccupationNumberProbabilityMeasurement.hpp":
    ctypedef vector[unsigned int] CppOccupationBounds "const std::vector<unsigned int>"

    cdef cppclass CppSubsystemOccupationNumberProbabilityMeasurement "SubsystemOccupationNumberProbabilityMeasurement" (CppBaseMeasurement):
        CppSubsystemOccupationNumberProbabilityMeasurement(unsigned int, shared_ptr[CppSubsystem]&)
        CppIntegerBinnedEstimate& get_estimate(CppOccupationBounds&)
        CppOccupationBounds& get_bounds()

cdef extern from "ParticleOperator.hpp":
    cdef cppclass CppSiteHop "SiteHop":
        CppSiteHop(CppLatticeSite& source, CppLatticeSite& destination, unsigned int species)

    bint is_valid_ParticleOperator "ParticleOperator::is_valid" (vector[CppSiteHop]&, CppLattice&, unsigned int N_species)

    cdef cppclass CppParticleOperator "ParticleOperator":
        CppParticleOperator(vector[CppSiteHop]&, shared_ptr[CppLattice]&)

cdef extern from "OperatorMeasurement.hpp":
    cdef cppclass CppOperatorMeasurement "OperatorMeasurement" (CppBaseMeasurement):
        CppOperatorMeasurement(unsigned int, CppParticleOperator&, CppBoundaryConditions&)
        CppComplexBinnedEstimate& get_estimate()


### implementation

from cython.operator cimport dereference as deref
from pyvmc.includes.libcpp.memory cimport auto_ptr

import collections

import numpy

from pyvmc.core cimport complex_t
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.boundary_conditions cimport BoundaryCondition
from pyvmc.core.lattice cimport Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.subsystem cimport Subsystem
from pyvmc.core.walk import StandardWalkPlan
from pyvmc.core.operator import SiteHop, BasicOperator, Operator
from pyvmc.core.estimate cimport Estimate_from_CppIntegerBinnedEstimate, Estimate_from_CppComplexBinnedEstimate
from pyvmc.core.measurement cimport BaseMeasurement
from pyvmc.core.measurement import BasicMeasurementPlan

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