from pyvmc.includes.boost.rational cimport rational as boost_rational

from pyvmc.core cimport lw_vector, complex_t

cdef extern from "BoundaryCondition.hpp":
    cdef cppclass CppBoundaryCondition "BoundaryCondition":
        CppBoundaryCondition()
        CppBoundaryCondition(boost_rational[int]&)

        bint is_initialized()

        bint operator==(const CppBoundaryCondition&)
        bint operator!=(const CppBoundaryCondition&)

        boost_rational[int] p()
        complex_t phase()

    CppBoundaryCondition cpp_periodic_bc "periodic_bc"
    CppBoundaryCondition cpp_antiperiodic_bc "antiperiodic_bc"
    CppBoundaryCondition cpp_open_bc "open_bc"

    ctypedef lw_vector[CppBoundaryCondition] CppBoundaryConditions "lw_vector<BoundaryCondition, MAX_DIMENSION>"

cdef class BoundaryCondition(object):
    cdef CppBoundaryCondition cpp
