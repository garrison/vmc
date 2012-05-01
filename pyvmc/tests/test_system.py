import json

import pytest

from pyvmc.system import LatticeSite, Lattice

def test_lattice_site():
    assert LatticeSite((1, 4)) == LatticeSite((1, 4), 0)

def test_lattice_iteration():
    lattice = Lattice((48, 4), 2)
    assert len(lattice) == (48 * 4 * 2)
    for i, s in enumerate(lattice):
        assert lattice[i] == s
        assert lattice.index(s) == i

    with pytest.raises(ValueError):
        lattice[len(lattice)]

def test_lattice_immutable():
    lattice = Lattice((48, 2))
    with pytest.raises(TypeError):
        lattice.x = 3
    with pytest.raises(TypeError):
        del lattice.x

def test_lattice_json():
    lattice = Lattice([24, 2])
    assert json.dumps(lattice.to_json()) == json.dumps({ "size": [24, 2] })
