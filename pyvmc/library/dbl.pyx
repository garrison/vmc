import numbers

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.orbitals import Orbitals
from pyvmc.core.orbitals cimport orbitals_to_orbitaldefinitions

class DBLWavefunction(Wavefunction):
    """$d$-wave Bose Liquid"""

    __slots__ = ('lattice', 'd1', 'd2', 'd1_exponent', 'd2_exponent')

    def init_validate(self, lattice, d1, d2, d1_exponent=1.0, d2_exponent=1.0):
        (lattice,) = super(DBLWavefunction, self).init_validate(lattice)
        assert isinstance(d1_exponent, numbers.Real)
        assert isinstance(d2_exponent, numbers.Real)
        return (
            lattice,
            Orbitals.from_description(d1, lattice),
            Orbitals.from_description(d2, lattice),
            float(d1_exponent),
            float(d2_exponent)
        )

    def to_json(self):
        return {
            'wavefunction': {
                'type': 'dbl',
                'orbitals-d1': self.d1.to_json(),
                'exponent-d1': self.d1_exponent,
                'orbitals-d2': self.d2.to_json(),
                'exponent-d2': self.d2_exponent,
            }
        }

    def to_wavefunction(self):
        cdef WavefunctionWrapper rv = WavefunctionWrapper()
        rv.sharedptr.reset(new CppDBLWavefunction(
            orbitals_to_orbitaldefinitions(self.d1, self.lattice),
            orbitals_to_orbitaldefinitions(self.d2, self.lattice),
            self.d1_exponent, self.d2_exponent
        ))
        return rv
