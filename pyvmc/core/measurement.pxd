from libcpp.vector cimport vector
from pyvmc.includes.boost.shared_ptr cimport shared_ptr

from pyvmc.core.lattice cimport CppLattice, CppLatticeSite
from pyvmc.core.boundary_conditions cimport CppBoundaryCondition, CppBoundaryConditions
from pyvmc.core.subsystem cimport CppSubsystem
from pyvmc.core.walk cimport CppWalk
from pyvmc.core cimport complex_t

cdef extern from "RunningEstimate.hpp":
    cdef cppclass CppIntegerRunningEstimate "RunningEstimate<unsigned int>":
        float get_recent_result()
        float get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

    cdef cppclass CppRealRunningEstimate "RunningEstimate<real_t>":
        float get_recent_result()
        float get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

    cdef cppclass CppComplexRunningEstimate "RunningEstimate<amplitude_t>":
        complex_t get_recent_result()
        complex_t get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

cdef extern from "BinnedEstimate.hpp":
    cdef cppclass CppIntegerBinnedEstimate "BinnedEstimate<unsigned int>" (CppIntegerRunningEstimate):
        pass

    cdef cppclass CppRealBinnedEstimate "BinnedEstimate<real_t>" (CppRealRunningEstimate):
        pass

    cdef cppclass CppComplexBinnedEstimate "BinnedEstimate<amplitude_t>" (CppComplexRunningEstimate):
        pass

cdef extern from "Measurement.hpp":
    cdef cppclass CppBaseMeasurement "BaseMeasurement":
        bint is_valid_walk(CppWalk&)

cdef class BaseMeasurement(object):
    cdef shared_ptr[CppBaseMeasurement] sharedptr

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
