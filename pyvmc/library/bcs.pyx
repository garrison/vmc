import numbers
import collections

from pyvmc.core.wavefunction import Wavefunction

class BCSWavefunction(Wavefunction):
    """Projected BCS wavefunction"""

    __slots__ = ('lattice', 'phi')

    def init_validate(self, lattice, phi):
        (lattice,) = super(BCSWavefunction, self).init_validate(lattice)

        assert isinstance(phi, collections.Sequence)
        phi = tuple(phi)
        assert len(phi) == len(lattice)
        assert all(isinstance(x, numbers.Complex) for x in phi)

        return (lattice, phi)

    @property
    def N_species(self):
        return 2

    def to_json(self):
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
            ('lattice', self.lattice.to_json()),
            ('phi', self.phi),
        ])

    def to_wavefunction(self):
        cdef vector[complex_t] phi
        for x in self.phi:
            phi.push_back(complex_t(x.real, x.imag))

        cdef WavefunctionWrapper rv = WavefunctionWrapper()
        rv.sharedptr.reset(new CppBCSWavefunction((<Lattice>self.lattice).sharedptr, phi))
        return rv
