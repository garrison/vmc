import numbers
import collections

from cython.operator cimport dereference as deref
from pyvmc.includes.libcpp.memory cimport unique_ptr

import numpy

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.utils import is_square_matrix

class ProjectedBCSWavefunction(Wavefunction):
    """Projected BCS wavefunction"""

    __slots__ = ('lattice', 'phi', 'N_up', 'N_dn')

    def init_validate(self, lattice, phi, N_up, N_dn=None):
        (lattice,) = super(ProjectedBCSWavefunction, self).init_validate(lattice)

        assert is_square_matrix(phi, len(lattice))
        if not all(numpy.isfinite(x) for x in numpy.nditer(phi)):
            raise RuntimeError("elements of phi matrix must be finite")

        assert isinstance(N_up, numbers.Integral)
        assert N_up > 0
        if N_dn is None:
            N_dn = N_up
        assert isinstance(N_dn, numbers.Integral)
        assert N_dn > 0
        assert N_up + N_dn <= len(lattice)

        return (lattice, phi, N_up, N_dn)

    def __hash__(self):
        # self.phi is a numpy.ndarray, which is unhashable, so we ignore it here
        return hash((self.lattice, self.N_up, self.N_dn))

    @property
    def N_species(self):
        return 2

    @property
    def N_filled(self):
        return (self.N_up, self.N_dn)

    def to_json(self):
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
            ('lattice', self.lattice.to_json()),
            #('phi', self.phi),  # phi is a huge matrix, so we store it in json_extra
            ('N_up', self.N_up),
            ('N_dn', self.N_dn),
        ])

    def to_json_extra(self):
        return collections.OrderedDict([
            ("phi", self.phi),
        ])

    def to_wavefunction(self):
        cdef unsigned int N_sites = len(self.lattice)
        phi = self.phi

        cdef unique_ptr[CppPhiMatrix] phimat
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
