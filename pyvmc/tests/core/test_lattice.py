import pytest

from pyvmc.utils import custom_json as json
from pyvmc.core.boundary_conditions import periodic, antiperiodic, open_bc
from pyvmc.core import LatticeSite, Lattice, HypercubicLattice, HexagonalLattice

def test_lattice_site():
    assert LatticeSite((1, 4)) == LatticeSite((1, 4), 0)

    assert LatticeSite((1, 4)) != LatticeSite((1, 4), 1)
    assert LatticeSite((0, 2)) != LatticeSite((0, 3))
    assert LatticeSite((0, 2)) != LatticeSite((0, 2, 3))

def test_lattice_site_inequality():
    lattice = Lattice([4, 4], 2)
    assert lattice[0] < lattice[1]
    assert not (lattice[5] < lattice[3])
    assert lattice[18] >= lattice[14]
    assert lattice[3] > lattice[1]
    assert lattice[7] <= lattice[19]

    assert lattice[7] <= lattice[7]
    assert lattice[5] >= lattice[5]
    assert not (lattice[3] < lattice[3])
    assert not (lattice[9] > lattice[9])

def test_lattice_site_displacement():
    assert LatticeSite([3, 1, -2]).displaced([-4, 1, 3]) == LatticeSite([-1, 2, 1])
    assert LatticeSite([0, 1], 1).negative_displaced([2, 3]) == LatticeSite([-2, -2], 1)
    with pytest.raises(ValueError):
        LatticeSite([3, 1, 2]).displaced([1, 2])
    with pytest.raises(AssertionError):
        LatticeSite([1, 3, 4]).displaced(LatticeSite([2, 3, 4]))

    assert LatticeSite([3, 1, -2]) + [-4, 1, 3] == LatticeSite([-1, 2, 1])
    assert LatticeSite([0, 1], 1) - [2, 3] == LatticeSite([-2, -2], 1)
    with pytest.raises(ValueError):
        LatticeSite([3, 1, 2]) + (1, 2)
    with pytest.raises(AssertionError):
        LatticeSite([1, 3, 4]) + LatticeSite([2, 3, 4])

    # addition must have the LatticeSite on the left
    with pytest.raises(TypeError):
        [1, 2] + LatticeSite([1, 3])
    with pytest.raises(TypeError):
        (1, 2) + LatticeSite([1, 3])

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

def test_lattice_bravais_site_iteration():
    lattice = Lattice((4, 8), 3)
    for site, bs in zip(lattice, lattice.iterate_bravais_sites()):
        assert site.bi == 0
        assert site.bs == bs
    assert len(list(lattice.iterate_bravais_sites())) == lattice.total_bravais_sites()

def test_lattice_immutable():
    lattice = Lattice((48, 2))
    with pytest.raises(Exception):
        lattice.x = 3
    with pytest.raises(Exception):
        del lattice.x
    with pytest.raises(Exception):
        lattice.dimensions[0] = 24

def test_lattice_json():
    lattice = Lattice([24, 2])
    assert json.dumps(lattice.to_json()) == '{"type": "Lattice", "dimensions": [24, 2], "basis_indices": 1}'
    assert Lattice.from_json(lattice.to_json()) == lattice

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
    assert lattice.enforce_boundary((24, 3), False) == (0, 3)
    assert lattice.enforce_boundary((-1, 4), False) == (23, 0)

    # LatticeSite
    assert lattice.enforce_boundary(LatticeSite((-1, 0)), False) == LatticeSite((23, 0))
    with pytest.raises(AssertionError):
        lattice.enforce_boundary(LatticeSite((-1, 0), 1), False)

def test_enforce_boundary_with_boundary_conditions():
    lattice = Lattice([8, 8], 2)

    # BravaisSite
    assert lattice.enforce_boundary((4, 4), (periodic, periodic)) == ((4, 4), 1)
    assert lattice.enforce_boundary((4, 4), (antiperiodic, antiperiodic)) == ((4, 4), 1)
    assert lattice.enforce_boundary((13, 4), (periodic, periodic)) == ((5, 4), 1)
    assert lattice.enforce_boundary((13, 4), (antiperiodic, periodic)) == ((5, 4), -1)
    assert lattice.enforce_boundary((-3, 4), (antiperiodic, periodic)) == ((5, 4), -1)
    assert lattice.enforce_boundary((-3, 4), (antiperiodic, open_bc)) == ((5, 4), -1)
    assert lattice.enforce_boundary((-3, 4), (open_bc, periodic)) == ((5, 4), 0)

    # LatticeSite
    assert lattice.enforce_boundary(LatticeSite((-3, 4), 1), (antiperiodic, periodic)) == (LatticeSite([5, 4], 1), -1)

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

class TestHexagonalSpace:
    def test_primitive_vectors(self):
        lattice = HexagonalLattice([8, 8])
        from numpy import pi, dot, abs
        for i, pv in enumerate(lattice.primitive_vectors):
            for j, rpv in enumerate(lattice.reciprocal_primitive_vectors):
                target = (2 * pi) if (i == j) else 0.0
                assert abs(dot(pv, rpv) - target) < .00000001

    def test_nearest_neighbors_distance(self):
        lattice = HexagonalLattice([8, 8])
        nn = lattice.nearest_neighbors(LatticeSite([0, 0]))
        nn = [lattice.spatial_coordinates(point) for point in nn]
        norms = [sum(a * a for a in point) for point in nn]
        assert min(norms) + .0000001 > max(norms)

    def test_second_nearest_neighbors_distance(self):
        lattice = HexagonalLattice([8, 8])
        nn = lattice.second_nearest_neighbors(LatticeSite([0, 0]))
        nn = [lattice.spatial_coordinates(point) for point in nn]
        norms = [sum(a * a for a in point) for point in nn]
        assert min(norms) + .0000001 > max(norms)

    def test_third_nearest_neighbors_distance(self):
        lattice = HexagonalLattice([8, 8])
        nn = lattice.third_nearest_neighbors(LatticeSite([0, 0]))
        nn = [lattice.spatial_coordinates(point) for point in nn]
        norms = [sum(a * a for a in point) for point in nn]
        assert min(norms) + .0000001 > max(norms)

    def test_basic_plaquettes_distance(self):
        lattice = HexagonalLattice([8, 8])
        for site1, site2, site3, site4 in lattice.basic_plaquettes():
            assert abs(sum(a * a for a in lattice.spatial_coordinates(site1 - site2.bs)) - 1) < .0000001
            assert abs(sum(a * a for a in lattice.spatial_coordinates(site2 - site3.bs)) - 1) < .0000001
            assert abs(sum(a * a for a in lattice.spatial_coordinates(site3 - site4.bs)) - 1) < .0000001
            assert abs(sum(a * a for a in lattice.spatial_coordinates(site4 - site1.bs)) - 1) < .0000001

class TestHypercubicSpace:
    def test_primitive_vectors(self):
        lattice = HypercubicLattice([8, 8, 16])
        from numpy import pi, dot, abs
        for i, pv in enumerate(lattice.primitive_vectors):
            for j, rpv in enumerate(lattice.reciprocal_primitive_vectors):
                target = (2 * pi) if (i == j) else 0.0
                assert abs(dot(pv, rpv) - target) < .00000001

    def test_nearest_neighbors_distance(self):
        lattice = HypercubicLattice([8, 16, 24])
        nn = lattice.nearest_neighbors(LatticeSite([0, 0, 0]))
        nn = [lattice.spatial_coordinates(point) for point in nn]
        norms = [sum(a * a for a in point) for point in nn]
        assert min(norms) + .0000001 > max(norms)
