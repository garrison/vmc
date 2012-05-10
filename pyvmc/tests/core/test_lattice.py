from fractions import Fraction

import pytest

from pyvmc.utils import complex_json as json
from pyvmc.core.boundary_conditions import periodic, antiperiodic
from pyvmc.core import LatticeSite, Lattice, HypercubicLattice, HexagonalLattice

def test_lattice_site():
    assert LatticeSite((1, 4)) == LatticeSite((1, 4), 0)

def test_lattice_constructor():
    with pytest.raises(AssertionError):
        Lattice([-1, 2])
    with pytest.raises(AssertionError):
        Lattice([])
    with pytest.raises(AssertionError):
        Lattice([0, 4])

def test_lattice_iteration():
    lattice = Lattice((48, 4), 2)
    assert len(lattice) == (48 * 4 * 2)
    for i, s in enumerate(lattice):
        assert lattice[i] == s
        assert lattice.index(s) == i
        assert s in lattice
    assert i + 1 == len(lattice)

    with pytest.raises(ValueError):
        lattice[len(lattice)]

    assert LatticeSite((48, 5)) not in lattice
    assert LatticeSite((46, 2), 2) not in lattice

    assert lattice.count(LatticeSite((5, 3))) == 1
    assert lattice.count(LatticeSite((245, 2))) == 0

def test_lattice_immutable():
    lattice = Lattice((48, 2))
    with pytest.raises(TypeError):
        lattice.x = 3
    with pytest.raises(TypeError):
        del lattice.x
    with pytest.raises(TypeError):
        lattice.dimensions[0] = 24

def test_lattice_json():
    lattice = Lattice([24, 2])
    assert json.dumps(lattice.to_json()) == json.dumps({ "size": [24, 2] })

def test_lattice_equality():
    lattice1 = Lattice([24, 2])
    lattice2 = Lattice([24, 2])
    lattice3 = Lattice([24, 2], 2)
    lattice4 = Lattice([48, 2])
    lattice5 = Lattice([48])

    assert lattice1 == lattice2
    assert not (lattice1 != lattice1)
    assert not (lattice1 != lattice2)
    assert lattice2 != lattice3
    assert lattice2 != lattice4
    assert lattice4 != lattice5

def test_lattice_repr():
    for lattice in Lattice([24, 2]), Lattice([8, 8], 2):
        assert lattice == eval(repr(lattice))

def test_vmc_core_lattice_correspondence():
    lattice = Lattice([32, 4], 2)
    assert lattice[3] == LatticeSite([3, 0], 0)
    assert lattice[32] == LatticeSite([0, 1], 0)

def test_enforce_boundary_without_boundary_conditions():
    lattice = Lattice([24, 4])

    # BravaisSite
    assert lattice.enforce_boundary((24, 3)) == (0, 3)
    assert lattice.enforce_boundary((-1, 4)) == (23, 0)

    # LatticeSite
    assert lattice.enforce_boundary(LatticeSite((-1, 0))) == LatticeSite((23, 0))
    with pytest.raises(AssertionError):
        lattice.enforce_boundary(LatticeSite((-1, 0), 1))

def test_enforce_boundary_with_boundary_conditions():
    lattice = Lattice([8, 8], 2)

    # BravaisSite
    assert lattice.enforce_boundary((4, 4), (periodic, periodic)) == ((4, 4), 0)
    assert lattice.enforce_boundary((4, 4), (antiperiodic, antiperiodic)) == ((4, 4), 0)
    assert lattice.enforce_boundary((13, 4), (periodic, periodic)) == ((5, 4), 0)
    assert lattice.enforce_boundary((13, 4), (antiperiodic, periodic)) == ((5, 4), Fraction(1, 2))
    assert lattice.enforce_boundary((-3, 4), (antiperiodic, periodic)) == ((5, 4), Fraction(1, 2))

    # LatticeSite
    assert lattice.enforce_boundary(LatticeSite((-3, 4), 1), (antiperiodic, periodic)) == (LatticeSite([5, 4], 1), Fraction(1, 2))

def test_hypercubic_lattice():
    lattice = HypercubicLattice([8, 8])
    assert eval(repr(lattice)) == lattice
    point = LatticeSite([0, 0])
    expected = set([
        LatticeSite([0, 1]),
        LatticeSite([1, 0]),
        LatticeSite([0, -1]),
        LatticeSite([-1, 0]),
    ])
    assert set(lattice.nearest_neighbors(point)) == expected
    assert len(lattice.nearest_neighbors(point)) == 2 * len(lattice.nearest_neighbors(point, False))

def test_hexagonal_lattice():
    lattice = HexagonalLattice([8, 8])
    assert eval(repr(lattice)) == lattice
    point = LatticeSite([4, 0])
    expected = set([
        LatticeSite([4, 1]),
        LatticeSite([5, 0]),
        LatticeSite([5, 1]),
        LatticeSite([4, -1]),
        LatticeSite([3, 0]),
        LatticeSite([3, -1]),
    ])
    assert set(lattice.nearest_neighbors(point)) == expected
    assert len(lattice.nearest_neighbors(point)) == 2 * len(lattice.nearest_neighbors(point, False))
