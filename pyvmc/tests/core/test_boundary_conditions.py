from fractions import Fraction

import pytest

from pyvmc.core.lattice import Lattice, LatticeSite
from pyvmc.core.boundary_conditions import periodic, antiperiodic, enforce_boundary

def test_boundary_condition_definitions():
    assert periodic == 1
    assert antiperiodic * 2 == 1

def test_enforce_boundary_without_boundary_conditions():
    lattice = Lattice([24, 4])

    # BravaisSite
    assert enforce_boundary((24, 3), lattice) == (0, 3)
    assert enforce_boundary((-1, 4), lattice) == (23, 0)

    # LatticeSite
    assert enforce_boundary(LatticeSite((-1, 0)), lattice) == LatticeSite((23, 0))
    with pytest.raises(AssertionError):
        enforce_boundary(LatticeSite((-1, 0), 1), lattice)

def test_enforce_boundary_with_boundary_conditions():
    lattice = Lattice([8, 8], 2)

    # BravaisSite
    assert enforce_boundary((4, 4), lattice, (periodic, periodic)) == ((4, 4), 0)
    assert enforce_boundary((4, 4), lattice, (antiperiodic, antiperiodic)) == ((4, 4), 0)
    assert enforce_boundary((13, 4), lattice, (periodic, periodic)) == ((5, 4), 0)
    assert enforce_boundary((13, 4), lattice, (antiperiodic, periodic)) == ((5, 4), Fraction(1, 2))
    assert enforce_boundary((-3, 4), lattice, (antiperiodic, periodic)) == ((5, 4), Fraction(1, 2))

    # LatticeSite
    assert enforce_boundary(LatticeSite((-3, 4), 1), lattice, (antiperiodic, periodic)) == (LatticeSite([5, 4], 1), Fraction(1, 2))
