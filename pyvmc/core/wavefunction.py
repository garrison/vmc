import abc

from pyvmc.core.lattice import Lattice
from pyvmc.core.orbitals import Orbitals
from pyvmc.utils.immutable import Immutable

class Wavefunction(Immutable):
    """Base class for all wavefunctions"""

    __slots__ = ('lattice',)

    def init_validate(self, lattice):
        assert isinstance(lattice, Lattice)
        return (lattice,)

    @abc.abstractmethod
    def to_json(self):
        return None

class FreeFermionWavefunction(Wavefunction):
    """Free fermion wavefunction, consists of a single determinant"""

    __slots__ = ('lattice', 'orbitals')

    def init_validate(self, lattice, orbitals):
        (lattice,) = super(FreeFermionWavefunction, self).init_validate(lattice)
        return lattice, Orbitals.from_description(orbitals, lattice)

    def to_json(self):
        return {
            'lattice': self.lattice.to_json(),
            'wavefunction': {
                'type': 'free-fermion',
                'orbitals': self.orbitals.to_json(),
            }
        }
