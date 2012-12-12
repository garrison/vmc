import numbers
import collections

from cython.operator cimport dereference as deref
from pyvmc.includes.libcpp.memory cimport auto_ptr

import numpy

from pyvmc.core.wavefunction import Wavefunction

class ProjectedBCSWavefunction(Wavefunction):
    """Projected BCS wavefunction"""

    __slots__ = ('lattice', 'phi', 'N_up', 'N_dn')

    def init_validate(self, lattice, phi, N_up, N_dn=None):
        (lattice,) = super(ProjectedBCSWavefunction, self).init_validate(lattice)

        assert isinstance(phi, numpy.ndarray)
        assert len(phi.shape) == 2
        assert phi.shape[0] == len(lattice)
        assert phi.shape[1] == len(lattice)

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
        assert False  # phi is a huge matrix now...
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
            ('lattice', self.lattice.to_json()),
            ('phi', self.phi),
            ('N_up', self.N_up),
            ('N_dn', self.N_dn),
        ])

    def to_wavefunction(self):
        cdef unsigned int N_sites = len(self.lattice)
        phi = self.phi

        cdef auto_ptr[CppPhiMatrix] phimat
        phimat.reset(new CppPhiMatrix(N_sites, N_sites))
        cdef unsigned int i, j
        for i in xrange(N_sites):
            for j in xrange(N_sites):
                c = phi[(i, j)]
                set_matrix_coeff(deref(phimat.get()), i, j, complex_t(c.real, c.imag))

        cdef WavefunctionWrapper rv = WavefunctionWrapper()
        assert self.N_up == self.N_dn
        rv.sharedptr.reset(new CppBCSWavefunction((<Lattice>self.lattice).sharedptr, deref(phimat), self.N_up))
        return rv