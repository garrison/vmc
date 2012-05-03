import json

import pytest

from pyvmc.core.lattice import Lattice
from pyvmc.core.orbitals import Orbitals, MomentaOrbitals, Bands
from pyvmc.core.boundary_conditions import periodic, antiperiodic

def test_momenta_orbitals():
    lattice = Lattice([8, 8])
    MomentaOrbitals(lattice, [(0, 0)], [periodic, periodic])
    with pytest.raises(AssertionError):
        MomentaOrbitals(lattice, [(9, 9)], [periodic, periodic])
    with pytest.raises(AssertionError):
        MomentaOrbitals(lattice, [(-1, 0)], [periodic, periodic])
    with pytest.raises(AssertionError):
        MomentaOrbitals(lattice, [(0, 0)], [periodic])
    with pytest.raises(AssertionError):
        MomentaOrbitals(lattice, [(0, 0, 0)], [periodic, periodic, periodic])
    with pytest.raises(AssertionError):
        MomentaOrbitals(lattice, [[0, 0]], [periodic, periodic])

def test_momenta_orbitals_to_json():
    lattice = Lattice([12, 2])
    momentum_sites = [(0, 0), (1, 0), (11, 0), (0, 1)]
    orbitals = MomentaOrbitals(lattice, momentum_sites, [periodic, periodic])
    assert json.dumps(orbitals.to_json()) == json.dumps({ 'filling': momentum_sites, 'boundary-conditions': [periodic, periodic] })

def test_single_band_orbitals():
    lattice = Lattice([8])
    assert set(Bands._single_band_orbitals((), 2, lattice)) == set([(0,), (7,)])

class TestBands:
    def test_single_band(self):
        lattice = Lattice([8])
        bands = Bands([4], [antiperiodic])
        assert Orbitals.from_description(bands, lattice) == MomentaOrbitals(lattice, [(0,), (1,), (7,), (6,)], [antiperiodic])

    def test_double_band(self):
        lattice = Lattice([8, 2])
        bands = Bands([3, 1], [periodic, periodic])
        assert Orbitals.from_description(bands, lattice) == MomentaOrbitals(lattice, [(0, 0), (1, 0), (7, 0), (0, 1)], [periodic, periodic])
