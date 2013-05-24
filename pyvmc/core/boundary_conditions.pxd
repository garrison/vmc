from pyvmc.includes.boost.rational cimport rational as boost_rational

from pyvmc.core cimport lw_vector, complex_t

cdef extern from "BoundaryCondition.hpp":
    cdef cppclass CppBoundaryCondition "BoundaryCondition<amplitude_t>":
        CppBoundaryCondition()
        CppBoundaryCondition(boost_rational[int]&)

        bint is_initialized()

        bint operator==(const CppBoundaryCondition&)
        bint operator!=(const CppBoundaryCondition&)

        boost_rational[int] p()
        complex_t phase()

    CppBoundaryCondition cpp_open_bc "BoundaryCondition<amplitude_t>::open"
    CppBoundaryCondition cpp_periodic_bc "BoundaryCondition<amplitude_t>::periodic"
    CppBoundaryCondition cpp_antiperiodic_bc "BoundaryCondition<amplitude_t>::antiperiodic"

    ctypedef lw_vector[CppBoundaryCondition] CppBoundaryConditions "lw_vector<BoundaryCondition<amplitude_t>, MAX_DIMENSION>"

cdef class BoundaryCondition(object):
    cdef CppBoundaryCondition cpp
