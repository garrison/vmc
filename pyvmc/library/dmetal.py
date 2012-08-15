import numbers

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.orbitals import Orbitals

class DMetalWavefunction(Wavefunction):
    """DMetal"""

    __slots__ = ("lattice", "d1", "d2", "f_up", "f_down", "d1_exponent", "d2_exponent", "f_up_exponent", "f_down_exponent")

    def init_validate(self, lattice, d1, d2, f_up, f_down, d1_exponent=1.0, d2_exponent=1.0, f_up_exponent=1.0, f_down_exponent=1.0):
        (lattice,) = super(DMetalWavefunction, self).init_validate(lattice)
        assert isinstance(d1_exponent, numbers.Real)
        assert isinstance(d2_exponent, numbers.Real)
        assert isinstance(f_up_exponent, numbers.Real)
        assert isinstance(f_down_exponent, numbers.Real)
        return (
            lattice,
            Orbitals.from_description(d1, lattice),
            Orbitals.from_description(d2, lattice),
            Orbitals.from_description(f_up, lattice),
            Orbitals.from_description(f_down, lattice),
            float(d1_exponent),
            float(d2_exponent),
            float(f_up_exponent),
            float(f_down_exponent)
        )

    def to_json(self):
        return {
            'wavefunction': {
                'type': 'dmetal',
                'orbitals-d1': self.d1.to_json(),
                'exponent-d1': self.d1_exponent,
                'orbitals-d2': self.d2.to_json(),
                'exponent-d2': self.d2_exponent,
                'orbitals-f_up': self.f_up.to_json(),
                'exponent-f_up': self.f_up_exponent,
                'orbitals-f_down': self.f_down.to_json(),
                'exponent-f_down': self.f_down_exponent,
            },
        }
