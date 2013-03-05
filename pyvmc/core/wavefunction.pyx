from cython.operator cimport dereference as deref
from pyvmc.includes.libcpp.memory cimport auto_ptr

import abc
import collections

import numpy

from pyvmc.core.lattice import Lattice
from pyvmc.core.orbitals import Orbitals
from pyvmc.core.orbitals cimport orbitals_to_orbitaldefinitions
from pyvmc.utils.immutable import Immutable
from pyvmc.utils import is_square_matrix, is_finite

class Wavefunction(Immutable):
    """Base class for all wavefunctions"""

    __slots__ = ('lattice',)

    def init_validate(self, lattice):
        assert isinstance(lattice, Lattice)
        return (lattice,)

    @abc.abstractmethod
    def to_json(self):
        return None

    def to_json_extra(self):
        # the values of this dict should be numpy arrays.  useful for saving
        # big things (such as the BCS phi matrix) that shouldn't go in json
        return {}

    @abc.abstractmethod
    def to_wavefunction(self):
        return None

    @abc.abstractproperty
    def N_species(self):
        raise NotImplementedError

    @abc.abstractproperty
    def N_filled(self):
        """should return a tuple with the filling of each species"""
        raise NotImplementedError

cdef shared_ptr[CppWavefunctionAmplitude] create_wfa(wf, RandomNumberGenerator rng) except *:
    assert rng.is_good()
    cdef Lattice lattice = wf.lattice
    cdef WavefunctionWrapper ww = wf.to_wavefunction()
    cdef shared_ptr[CppWavefunctionAmplitude] wfa = ww.sharedptr.get().create_nonzero_wavefunctionamplitude(ww.sharedptr, deref(rng.autoptr.get()))
    if wfa.get() is NULL:
        raise RuntimeError("could not find a nonzero wavefunction configuration")
    return wfa

class FreeFermionWavefunction(Wavefunction):
    """Free fermion wavefunction, consists of a single determinant"""

    __slots__ = ('lattice', 'orbitals', 'jastrow')

    def init_validate(self, lattice, orbitals, jastrow=None):
        (lattice,) = super(FreeFermionWavefunction, self).init_validate(lattice)
        assert isinstance(orbitals, collections.Sequence)
        assert len(orbitals) > 0
        orbitals = tuple([Orbitals.from_description(orb, lattice) for orb in orbitals])
        if jastrow is not None:
            assert isinstance(jastrow, JastrowFactor)
        return lattice, orbitals, jastrow

    @property
    def N_species(self):
        return len(self.orbitals)

    @property
    def N_filled(self):
        return tuple([len(orbs) for orbs in self.orbitals])

    def to_json(self):
        d = [
            ('type', self.__class__.__name__),
            ('orbitals', [orbitals.to_json() for orbitals in self.orbitals]),
        ]
        if self.jastrow is not None:
            d.append(('jastrow', self.jastrow.to_json()))
        return collections.OrderedDict(d)

    def to_wavefunction(self):
        cdef vector[shared_ptr[const_CppOrbitalDefinitions]] orbital_defs
        for orbitals in self.orbitals:
            orbital_defs.push_back(orbitals_to_orbitaldefinitions(orbitals, self.lattice))

        cdef shared_ptr[CppJastrowFactor] jastrow_sharedptr
        if self.jastrow is not None:
            jastrow_sharedptr = (<JastrowFactor>self.jastrow).sharedptr

        cdef WavefunctionWrapper rv = WavefunctionWrapper()
        rv.sharedptr.reset(new CppFreeFermionWavefunction(orbital_defs, jastrow_sharedptr))
        return rv

cdef class NoDoubleOccupancyProjector(JastrowFactor):
    def __init__(self):
        self.sharedptr.reset(new CppNoDoubleOccupancyProjector())

    def to_json(self):
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
        ])

    def __richcmp__(self, other, int op):
        if op == 2:  # ==
            return (self.__class__ == other.__class__)
        elif op == 3:  # !=
            return (self.__class__ != other.__class__)
        # we don't implement <, <=, >, >=
        raise NotImplementedError

    def __hash__(self):
        return hash(self.__class__.__name__)

    def __repr__(self):
        return "{}()".format(self.__class__.__name__)

# for compatibility with previous versions
SingleOccupancyProjector = NoDoubleOccupancyProjector

cdef class TwoBodyJastrowFactor(JastrowFactor):
    """$\exp ( -\frac{1}{2} \sum_{ij} u_{ij} n_i n_j )$

    where $n_i$ counts the number of particles of all species on a site
    """

    cdef object _correlation_matrix

    def __init__(self, correlation_matrix):
        assert is_square_matrix(correlation_matrix)
        N_sites = correlation_matrix.shape[0]  # since we don't have access to the lattice here
        if not all(is_finite(x) for x in numpy.nditer(correlation_matrix)):
            raise RuntimeError("elements of correlation matrix must be finite")

        self._correlation_matrix = correlation_matrix

        # CYTHON-LIMITATION: we are forced to use an auto_ptr here because
        # Cython doesn't allow us to pass arguments to the constructor
        cdef auto_ptr[CppTwoBodyJastrowMatrix] corrmat
        corrmat.reset(new CppTwoBodyJastrowMatrix(N_sites, N_sites))
        cdef unsigned int i, j
        for i in xrange(N_sites):
            for j in xrange(N_sites):
                set_matrix_coeff(deref(corrmat.get()), i, j, correlation_matrix[(i, j)])

        self.sharedptr.reset(new CppTwoBodyJastrowFactor(deref(corrmat)))

    def to_json(self):
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
            ('correlation_matrix', [list(row) for row in self._correlation_matrix]),
        ])

    def __richcmp__(self, other, int op):
        if op == 2:  # ==
            return (self.__class__ == other.__class__ and
                    self._correlation_matrix == other._correlation_matrix)
        elif op == 3:  # !=
            return (self.__class__ != other.__class__ or
                    self._correlation_matrix != other._correlation_matrix)
        # we don't implement <, <=, >, >=
        raise NotImplementedError

    def __hash__(self):
        # FIXME: this is not very effective... see
        # http://stackoverflow.com/questions/5386694/fast-way-to-hash-numpy-objects-for-caching
        # for a way to also hash on self._correlation_matrix
        return hash(self.__class__.__name__)
