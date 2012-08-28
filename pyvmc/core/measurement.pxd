from libcpp.vector cimport vector

from pyvmc.boost.shared_ptr cimport shared_ptr
from pyvmc.boost.rational cimport rational as boost_rational
from pyvmc.core.lattice cimport CppLattice, CppLatticeSite, lw_vector
from pyvmc.core.subsystem cimport CppSubsystem

cdef extern from "Measurement.hpp":
    cdef cppclass CppBaseMeasurement "BaseMeasurement":
        pass

cdef class BaseMeasurement(object):
    cdef shared_ptr[CppBaseMeasurement] *sharedptr

cdef extern from "SubsystemOccupationNumberProbabilityMeasurement.hpp":
    cdef cppclass CppSubsystemOccupationNumberProbabilityMeasurement "SubsystemOccupationNumberProbabilityMeasurement" (CppBaseMeasurement):
        CppSubsystemOccupationNumberProbabilityMeasurement(int, shared_ptr[CppSubsystem])

cdef extern from "BoundaryCondition.hpp":
    cdef cppclass CppBoundaryCondition "BoundaryCondition":
        CppBoundaryCondition(boost_rational[int]&)

    CppBoundaryCondition periodic_bc
    CppBoundaryCondition antiperiodic_bc

    ctypedef lw_vector[CppBoundaryCondition] CppBoundaryConditions "lw_vector<BoundaryCondition, MAX_DIMENSION>"

cdef extern from "ParticleOperator.hpp":
    cdef cppclass CppSiteHop "SiteHop":
        CppSiteHop(CppLatticeSite& source, CppLatticeSite& destination, int species)

    bint is_valid_ParticleOperator "ParticleOperator::is_valid" (vector[CppSiteHop]&, CppLattice&, int N_species)

    cdef cppclass CppParticleOperator "ParticleOperator":
        CppParticleOperator(vector[CppSiteHop]&, shared_ptr[CppLattice]&)

cdef extern from "OperatorMeasurement.hpp":
    cdef cppclass CppOperatorMeasurement "OperatorMeasurement" (CppBaseMeasurement):
        CppOperatorMeasurement(int, CppParticleOperator&, bint, CppBoundaryConditions*)
