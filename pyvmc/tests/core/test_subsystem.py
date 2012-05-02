import json

import pytest

from pyvmc.core import Lattice, SimpleSubsystem

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

def test_subsystem_equality():
    lattice1 = Lattice([24, 4])
    lattice2 = Lattice([24, 4])
    lattice3 = Lattice([16, 4])

    subsystem1 = SimpleSubsystem([4, 4], lattice1)
    subsystem2 = SimpleSubsystem([4, 4], lattice2)
    subsystem3 = SimpleSubsystem([4, 4], lattice3)
    subsystem4 = SimpleSubsystem([3, 4], lattice1)

    assert subsystem1 == subsystem2
    assert not (subsystem1 != subsystem2)
    assert subsystem2 != subsystem3
    assert not (subsystem2 == subsystem3)
    assert subsystem1 != subsystem4

def test_subsystem_repr():
    lattice = Lattice([24, 4])
    for subsystem in SimpleSubsystem([4, 4], lattice), SimpleSubsystem([8, 4], lattice):
        assert subsystem == eval(repr(subsystem))
