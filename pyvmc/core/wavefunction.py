import abc

from pyvmc.core.lattice import Lattice
from pyvmc.core.orbitals import Orbitals

class Wavefunction(object):
    """Base class for all wavefunctions"""

    __metaclass__ = abc.ABCMeta

    def __init__(self, lattice):
        assert isinstance(lattice, Lattice)
        self.lattice = lattice

    @abc.abstractmethod
    def to_json(self):
        return None

class FreeFermionWavefunction(Wavefunction):
    """Free fermion wavefunction, consists of a single determinant"""

    def __init__(self, lattice, orbitals):
        super(FreeFermionWavefunction, self).__init__(lattice)
        object.__setattr__(self, "orbitals", Orbitals.from_description(orbitals, lattice))

    def to_json(self):
        return {
            'lattice': self.lattice.to_json(),
            'wavefunction': {
                'type': 'free-fermion',
                'orbitals': self.orbitals.to_json(),
            }
        }
