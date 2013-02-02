from fractions import Fraction
import collections
import numbers

periodic_bc = periodic = 1
antiperiodic_bc = antiperiodic = Fraction(1, 2)
open_bc = 0

def boundary_condition_to_string(bc):
    return {
        periodic: 'periodic',
        antiperiodic: 'antiperiodic',
        open_bc: 'open',
    }.get(bc, repr(bc))

def valid_boundary_conditions(boundary_conditions, n_dimensions):
    assert isinstance(n_dimensions, numbers.Integral) and n_dimensions > 0
    return bool(isinstance(boundary_conditions, collections.Sequence) and
                len(boundary_conditions) == n_dimensions and
                all(isinstance(bc, numbers.Real) and 0 <= bc <= 1
                    for bc in boundary_conditions))

def valid_closed_boundary_conditions(boundary_conditions, n_dimensions):
    return bool(valid_boundary_conditions(boundary_conditions, n_dimensions) and
                all(bc != 0 for bc in boundary_conditions))
