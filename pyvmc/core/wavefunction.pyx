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

cdef shared_ptr[CppWavefunctionAmplitude] create_wfa(wf):
    from pyvmc.utils import complex_json as json
    cdef unicode input_unicode = unicode(json.dumps(wf.to_json()))
    cdef bytes input_bytes = input_unicode.encode('UTF-8')
    cdef char* input_cstr = input_bytes
    cdef Lattice lattice = wf.lattice
    rng = RandomNumberGenerator()
    return create_wfa_from_json(input_cstr, lattice.sharedptr, rng.autoptr)

class FreeFermionWavefunction(Wavefunction):
    """Free fermion wavefunction, consists of a single determinant"""

    __slots__ = ('lattice', 'orbitals')

    def init_validate(self, lattice, orbitals):
        (lattice,) = super(FreeFermionWavefunction, self).init_validate(lattice)
        return lattice, Orbitals.from_description(orbitals, lattice)

    def to_json(self):
        return {
            'wavefunction': {
                'type': 'free-fermion',
                'orbitals': self.orbitals.to_json(),
            }
        }
