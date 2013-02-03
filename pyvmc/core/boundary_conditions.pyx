from fractions import Fraction
import collections
import numbers

from pyvmc.includes.boost.rational cimport rational as boost_rational
from pyvmc.core.lattice cimport MAX_DIMENSION

cdef class BoundaryCondition(object):
    """represents a boundary condition in a single direction

    (e.g. periodic, antiperiodic, open, twisted)

    remains constant after initialization
    """

    def __init__(self, bc):
        assert isinstance(bc, numbers.Rational)
        assert 0 <= bc <= 1
        self.cpp = CppBoundaryCondition(boost_rational[int](bc.numerator, bc.denominator))

    property p:
        def __get__(self):
            return Fraction(self.cpp.p().numerator(), self.cpp.p().denominator())

    property phase:
        def __get__(self):
            return complex(self.cpp.phase().real(), self.cpp.phase().imag())

    def __richcmp__(BoundaryCondition self not None, other, int op):
        if self.__class__ != other.__class__:
            if op == 2:  #  ==
                return False
            elif op == 3: #  !=
                return True
            else:
                raise NotImplementedError

        cdef BoundaryCondition other_ = other

        if op == 2:  # ==
            return self.cpp == other_.cpp
        elif op == 3:  # !=
            return self.cpp != other_.cpp
        else:
            raise NotImplementedError

    def __hash__(self):
        if not self.cpp.is_initialized():
            # FIXME: why is it even possible for this to be uninitialized?
            return 0
        return 256 * self.cpp.p().numerator() + self.cpp.p().denominator()

    def __repr__(self):
        if self == periodic_bc:
            return 'periodic_bc'
        elif self == antiperiodic_bc:
            return 'antiperiodic_bc'
        elif self == open_bc:
            return 'open_bc'
        else:
            return "{}({})".format(self.__class__.__name__, repr(self.p))

    def __str__(self):
        if self == periodic_bc:
            return 'periodic'
        elif self == antiperiodic_bc:
            return 'antiperiodic'
        elif self == open_bc:
            return 'open'
        else:
            return repr(self.p)

collections.Hashable.register(BoundaryCondition)

cdef BoundaryCondition_from_cpp(CppBoundaryCondition cppbc):
    assert cppbc.is_initialized()
    cdef BoundaryCondition bc = BoundaryCondition.__new__(BoundaryCondition)
    bc.cpp = cppbc
    return bc

periodic_bc = periodic = BoundaryCondition_from_cpp(cpp_periodic_bc)
antiperiodic_bc = antiperiodic = BoundaryCondition_from_cpp(cpp_antiperiodic_bc)
open_bc = BoundaryCondition_from_cpp(cpp_open_bc)

def boundary_condition_to_string(bc):
    from warnings import warn
    warn("boundary_condition_to_string() is deprecated.  use str(bc) instead.", DeprecationWarning, stacklevel=2)
    return str(bc)

def valid_boundary_conditions(boundary_conditions, n_dimensions):
    assert isinstance(n_dimensions, numbers.Integral) and 0 < n_dimensions <= MAX_DIMENSION
    return bool(isinstance(boundary_conditions, collections.Sequence) and
                len(boundary_conditions) == n_dimensions and
                all(isinstance(bc, BoundaryCondition)
                    for bc in boundary_conditions))

def valid_closed_boundary_conditions(boundary_conditions, n_dimensions):
    return bool(valid_boundary_conditions(boundary_conditions, n_dimensions) and
                all(bc != open_bc for bc in boundary_conditions))
