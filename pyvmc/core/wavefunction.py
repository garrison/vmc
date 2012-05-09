import abc
import collections

from pyvmc.core.lattice import Lattice
from pyvmc.core.orbitals import Orbitals

class Wavefunction(collections.Hashable):
    """Base class for all wavefunctions"""

    __metaclass__ = abc.ABCMeta

    __slots__ = ('lattice',)

    def __init__(self, lattice):
        assert isinstance(lattice, Lattice)
        object.__setattr__(self, "lattice", lattice)

    @abc.abstractmethod
    def to_json(self):
        return None

    def __setattr__(self, name, value):
        raise TypeError

    def __delattr__(self, name):
        raise TypeError

class FreeFermionWavefunction(Wavefunction):
    """Free fermion wavefunction, consists of a single determinant"""

    __slots__ = ('lattice', 'orbitals')

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

    def __hash__(self):
        return hash(self.lattice) | hash(self.orbitals)
