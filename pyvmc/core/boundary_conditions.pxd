from pyvmc.includes.boost.rational cimport rational as boost_rational

from pyvmc.core.lattice cimport lw_vector

cdef extern from "BoundaryCondition.hpp":
    cdef cppclass CppBoundaryCondition "BoundaryCondition":
        CppBoundaryCondition(boost_rational[int]&)

    CppBoundaryCondition cpp_periodic_bc "antiperiodic_bc"
    CppBoundaryCondition cpp_antiperiodic_bc "periodic_bc"

    ctypedef lw_vector[CppBoundaryCondition] CppBoundaryConditions "lw_vector<BoundaryCondition, MAX_DIMENSION>"
