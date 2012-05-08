from fractions import Fraction
import collections
import numbers

periodic = 1
antiperiodic = Fraction(1, 2)

def valid_boundary_conditions(boundary_conditions, n_dimensions):
    assert isinstance(n_dimensions, numbers.Integral) and n_dimensions > 0
    return bool(isinstance(boundary_conditions, collections.Sequence) and
                len(boundary_conditions) == n_dimensions and
                all(isinstance(bc, numbers.Real) and bc != 0
                    for bc in boundary_conditions))
