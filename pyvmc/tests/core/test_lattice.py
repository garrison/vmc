import json

import pytest

from pyvmc.core import LatticeSite, Lattice

def test_lattice_site():
    assert LatticeSite((1, 4)) == LatticeSite((1, 4), 0)

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

def test_lattice_immutable():
    lattice = Lattice((48, 2))
    with pytest.raises(TypeError):
        lattice.x = 3
    with pytest.raises(TypeError):
        del lattice.x

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

def test_vmc_core_lattice_correspondence():
    lattice = Lattice([32, 4], 2)
    assert lattice[3] == LatticeSite([3, 0], 0)
    assert lattice[32] == LatticeSite([0, 1], 0)
