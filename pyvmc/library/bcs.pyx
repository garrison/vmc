import numbers
import collections

from pyvmc.core.wavefunction import Wavefunction

class ProjectedBCSWavefunction(Wavefunction):
    """Projected BCS wavefunction"""

    __slots__ = ('lattice', 'phi', 'N_up', 'N_dn')

    def init_validate(self, lattice, phi, N_up, N_dn=None):
        (lattice,) = super(ProjectedBCSWavefunction, self).init_validate(lattice)

        assert isinstance(phi, collections.Sequence)
        phi = tuple(phi)
        assert len(phi) == len(lattice)
        assert all(isinstance(x, numbers.Complex) for x in phi)

        assert isinstance(N_up, numbers.Integral)
        assert N_up > 0
        if N_dn is None:
            N_dn = N_up
        assert isinstance(N_dn, numbers.Integral)
        assert N_dn > 0
        assert N_up + N_dn <= len(lattice)

        return (lattice, phi, N_up, N_dn)

    @property
    def N_species(self):
        return 2

    def to_json(self):
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
            ('lattice', self.lattice.to_json()),
            ('phi', self.phi),
            ('N_up', self.N_up),
            ('N_dn', self.N_dn),
        ])

    def to_wavefunction(self):
        cdef vector[complex_t] phi
        for x in self.phi:
            phi.push_back(complex_t(x.real, x.imag))

        cdef WavefunctionWrapper rv = WavefunctionWrapper()
        assert self.N_up == self.N_dn
        rv.sharedptr.reset(new CppBCSWavefunction((<Lattice>self.lattice).sharedptr, phi, self.N_up))
        return rv
