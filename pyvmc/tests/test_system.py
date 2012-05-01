import json

import pytest

from pyvmc.system import LatticeSite, Lattice, SimpleSubsystem

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

def test_vmc_core_lattice_correspondence():
    lattice = Lattice([32, 4], 2)
    assert lattice[3] == LatticeSite([3, 0], 0)
    assert lattice[32] == LatticeSite([0, 1], 0)

def test_simple_subsystem():
    lattice = Lattice([8, 8], 2)
    subsystem = SimpleSubsystem([4, 6], lattice)

    inside = set()
    for site in lattice:
        if site in subsystem:
            inside.add(site)
        else:
            with pytest.raises(ValueError):
                subsystem.index(site)

    assert len(inside) == len(subsystem)

    for i, site in enumerate(subsystem):
        assert subsystem[i] == site
        assert subsystem.index(site) == i
        inside.remove(site)

    assert len(inside) == 0

def test_simple_subsystem_json():
    lattice = Lattice([16, 8], 2)
    subsystem = SimpleSubsystem([4, 4], lattice)
    assert json.dumps(subsystem.to_json()) == json.dumps({ 'type': 'simple', 'dimensions': (4, 4) })
