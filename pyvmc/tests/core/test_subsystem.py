import pytest

from pyvmc.utils import custom_json as json
from pyvmc.core.lattice import Lattice
from pyvmc.core.subsystem import Subsystem, SimpleSubsystem, CustomSubsystem

class TestSimpleSubsystem:
    def test_simple_subsystem(self):
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

    def test_json(self):
        lattice = Lattice([16, 8], 2)
        subsystem = SimpleSubsystem([4, 4], lattice)
        assert json.dumps(subsystem.to_json()) == '{"type": "SimpleSubsystem", "dimensions": [4, 4]}'

        assert Subsystem.from_json(subsystem.to_json(), lattice) == subsystem
        assert SimpleSubsystem.from_json(subsystem.to_json(), lattice) == subsystem

    def test_equality(self):
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

    def test_repr(self):
        lattice = Lattice([24, 4])
        for subsystem in SimpleSubsystem([4, 4], lattice), SimpleSubsystem([8, 4], lattice):
            assert subsystem == eval(repr(subsystem))

class TestCustomSubsystem:
    def test_custom_subsystem(self):
        lattice = Lattice([8, 8], 2)
        subsystem = CustomSubsystem(lambda site: site.bs[0] % 2 == 0, lattice)

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

    def test_json(self):
        lattice = Lattice([16, 8], 2)
        subsystem = CustomSubsystem(lambda site: site.bs[0] == 0 and site.bs[1] < 4, lattice)
        assert json.dumps(subsystem.to_json()) == '{"type": "CustomSubsystem", "_site_indices": [0, 16, 32, 48, 128, 144, 160, 176]}'

        assert Subsystem.from_json(subsystem.to_json(), lattice) == subsystem
        assert CustomSubsystem.from_json(subsystem.to_json(), lattice) == subsystem

    def test_equality(self):
        lattice1 = Lattice([24, 4])
        lattice2 = Lattice([24, 4])
        lattice3 = Lattice([16, 4])

        subsystem1 = CustomSubsystem(lambda site: site.bs[0] == 0, lattice1)
        subsystem2 = CustomSubsystem(lambda site: site.bs[0] == 0, lattice2)
        subsystem3 = CustomSubsystem(lambda site: site.bs[0] == 0, lattice3)
        subsystem4 = CustomSubsystem(lambda site: site.bs[0] < 2, lattice1)

        assert subsystem1 == subsystem2
        assert not (subsystem1 != subsystem2)
        assert subsystem2 != subsystem3
        assert not (subsystem2 == subsystem3)
        assert subsystem1 != subsystem4

    def test_repr(self):
        lattice = Lattice([24, 4])
        for subsystem in CustomSubsystem(lambda site: site.bs == (4, 3), lattice), CustomSubsystem(lambda site: site.bs[0] % 2 == 0, lattice):
            assert subsystem == eval(repr(subsystem))
