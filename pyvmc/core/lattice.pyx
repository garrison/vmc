"""
"""

import abc
import numbers
import collections

import numpy
from cython.operator cimport dereference as deref

from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.utils import product
from pyvmc.utils.immutable import Immutable
from pyvmc.constants import two_pi, sqrt_three_over_two, pi

cdef extern from "Lattice.hpp":
    cdef unsigned int MAX_DIMENSION

    cdef cppclass CppLatticeSite "LatticeSite":
        CppLatticeSite(unsigned int)
        int operator[](int)
        void set_bs_coordinate(int, int)
        int n_dimensions()
        int basis_index

cdef class LatticeSite(object):
    """represents a site on a lattice

    bs == bravais site
    bi == basis index

    remains constant after initialization
    """

    cdef CppLatticeSite *thisptr

    def __cinit__(self, bs, bi=0):
        assert isinstance(bs, collections.Sequence)
        assert len(bs) > 0
        assert all(isinstance(x, numbers.Integral) for x in bs)
        assert isinstance(bi, numbers.Integral) and bi >= 0
        if len(bs) > MAX_DIMENSION:
            raise ValueError("provided site has greater than {} dimensions".format(MAX_DIMENSION))
        self.thisptr = new CppLatticeSite(len(bs))
        cdef int i, x
        for i, x in enumerate(bs):
            self.thisptr.set_bs_coordinate(i, x)
        self.thisptr.basis_index = bi

    property bs:
        def __get__(self):
            cdef int i
            return tuple([deref(self.thisptr)[i] for i in range(self.thisptr.n_dimensions())])

    property bi:
        def __get__(self):
            return self.thisptr.basis_index

    def __hash__(self):
        return hash(self.bs) | hash(self.bi)

    def __richcmp__(self, other, int op):
        if op == 2:  # ==
            return (self.__class__ == other.__class__ and
                    self.bs == other.bs and
                    self.bi == other.bi)
        elif op == 3:  # !=
            return (self.__class__ != other.__class__ or
                    self.bs != other.bs or
                    self.bi != other.bi)
        # we don't implement <, <=, >, >=
        raise NotImplementedError

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.bs),
                                   self.bi)

collections.Hashable.register(LatticeSite)

class Lattice(Immutable):
    __slots__ = ('dimensions', 'basis_indices')

    def init_validate(self, dimensions, basis_indices=1):
        assert isinstance(dimensions, collections.Sequence)
        assert len(dimensions) != 0
        dimensions = tuple(dimensions)
        assert all(isinstance(x, numbers.Integral) and x > 0 for x in dimensions)
        assert isinstance(basis_indices, numbers.Integral) and basis_indices > 0
        return (dimensions, basis_indices)

    def to_json(self):
        assert self.basis_indices == 1  # for now
        return {
            'size': self.dimensions,
        }

    def abstract_lattice(self):
        """This will always return a Lattice object, never a LatticeRealization"""
        return self

    def enforce_boundary(self, site, boundary_conditions=None):
        """Enforce the boundary of a site which may be outside the lattice

        Any quantum amplitudes should be multiplied by
        exp(2 * pi * i * phase_adjustment)
        """
        lattice_dimensions = self.dimensions
        is_LatticeSite = isinstance(site, LatticeSite)
        bravais_site = site.bs if is_LatticeSite else site
        assert isinstance(bravais_site, collections.Sequence)
        assert len(bravais_site) == len(lattice_dimensions)
        assert all(isinstance(x, numbers.Integral) for x in bravais_site)
        new_bravais_site = tuple(x % length for x, length in zip(bravais_site, lattice_dimensions))
        new_site = LatticeSite(new_bravais_site, site.bi) if is_LatticeSite else new_bravais_site
        assert (not is_LatticeSite) or new_site in self
        if boundary_conditions is None:
            return new_site
        else:
            assert valid_boundary_conditions(boundary_conditions, len(lattice_dimensions))
            phase_adjustment = sum((x // length) * bc for x, length, bc
                                   in zip(bravais_site, lattice_dimensions, boundary_conditions)) % 1
            return new_site, phase_adjustment

    def __len__(self):
        return product(self.dimensions) * self.basis_indices

    def __iter__(self):
        dimensions = self.dimensions
        d = len(dimensions)
        s_bound = dimensions + (self.basis_indices,)
        s_len = d + 1
        s = [0] * s_len
        while True:
            yield LatticeSite(s[:d], s[-1])
            for i in xrange(s_len):
                s[i] += 1
                if s[i] < s_bound[i]:
                    break
                s[i] = 0
            else:
                raise StopIteration

    def __getitem__(self, index):
        assert isinstance(index, numbers.Integral)
        if index < 0 or index >= len(self):
            raise ValueError
        dimensions = self.dimensions
        bs = []
        for d in dimensions:
            bs.append(index % d)
            index //= d
        # "index" now represents the bravais index
        site = LatticeSite(bs, index)
        assert site in self
        return site

    def index(self, site):
        if site not in self:
            raise ValueError
        dimensions = self.dimensions
        index = 0
        offset = 1
        for i, x in enumerate(site.bs):
            index += x * offset
            offset *= dimensions[i]
        index += site.bi * offset
        assert index < len(self)
        return index

    def __contains__(self, site):
        assert isinstance(site, LatticeSite)
        return bool(len(site.bs) == len(self.dimensions) and
                    all(x < y for x, y in zip(site.bs, self.dimensions)) and
                    site.bi < self.basis_indices)

    def count(self, x):
        return 1 if x in self else 0

collections.Hashable.register(Lattice)
collections.Sequence.register(Lattice)

class LatticeRealization(Lattice):
    __metaclass__ = abc.ABCMeta

    __slots__ = ('dimensions', 'basis_indices')

    def abstract_lattice(self):
        return Lattice(self.dimensions, self.basis_indices)

    @abc.abstractproperty
    def primitive_vectors(self):
        return None

    @abc.abstractproperty
    def reciprocal_primitive_vectors(self):
        return None

    @property
    def basis_offsets(self):
        assert self.basis_indices == 1
        return (
            tuple(0.0 for i in xrange(len(self.dimensions))),
        )

    def real_space_point(self, point):
        # NOTE: this method doesn't even require that the resulting point is on
        # the lattice
        assert isinstance(point, LatticeSite)
        basis_offsets = self.basis_offsets
        assert len(basis_offsets) == self.basis_indices
        assert point.bi < self.basis_indices
        return tuple(sum(z) for z in zip(basis_offsets[point.bi],
                                         *(numpy.multiply(p, pv) for p, pv in zip(point.bs, self.primitive_vectors))))

    @abc.abstractmethod
    def nearest_neighbors(self, point, double_count=True):
        """Returns a point's nearest neighbors.

        If double_count evaluates to True (the default), all nearest neighbors
        are returned.  If double_count is False, this method will
        systematically return only half of the nearest neighbors for the point,
        such that summing over all points' neighest neighbors will sum over
        each bond precisely once.

        This method will, in some instances, return sites that aren't even on
        the lattice.  It is up to the user of this method to choose whether to
        call enforce_boundary() (e.g. in the case of periodic or twisted
        boundary conditions) or to throw such sites away (e.g. in the case of
        open boundary conditions).
        """
        assert point in self
        raise NotImplementedError

    def __repr__(self):
        return "%s(%s)" % (
            self.__class__.__name__,
            repr(self.dimensions)
        )

class HypercubicLattice(LatticeRealization):
    __slots__ = ('dimensions', 'basis_indices')

    a = 1.0   # lattice spacing

    def init_validate(self, dimensions):
        return super(HypercubicLattice, self).init_validate(dimensions, 1)

    @property
    def primitive_vectors(self):
        rv = []
        num_dimensions = len(self.dimensions)
        for d in xrange(num_dimensions):
            v = [0.0] * num_dimensions
            v[d] = self.a
            rv.append(tuple(v))
        return rv

    @property
    def reciprocal_primitive_vectors(self):
        rv = []
        b = two_pi / self.a
        num_dimensions = len(self.dimensions)
        for d in xrange(num_dimensions):
            v = [0.0] * num_dimensions
            v[d] = b
            rv.append(tuple(v))
        return rv

    def nearest_neighbors(self, point, double_count=True):
        assert point in self
        bs = point.bs
        rv = []
        for d in xrange(len(self.dimensions)):
            p = list(bs)
            p[d] += 1
            rv.append(LatticeSite(p, 0))
            if double_count:
                p[d] -= 2
                rv.append(LatticeSite(p, 0))
        return rv

class HexagonalLattice(LatticeRealization):
    __slots__ = ('dimensions', 'basis_indices')

    a = 1.0   # lattice spacing

    def init_validate(self, dimensions):
        assert len(dimensions) == 2
        return super(HexagonalLattice, self).init_validate(dimensions, 1)

    @property
    def primitive_vectors(self):
        a = self.a
        return (
            (a, 0.0),
            (-.5 * a, sqrt_three_over_two * a),
        )

    @property
    def reciprocal_primitive_vectors(self):
        a = self.a
        return (
            (two_pi / self.a, pi / sqrt_three_over_two / self.a),
            (0.0, two_pi / sqrt_three_over_two / self.a),
        )

    def nearest_neighbors(self, point, double_count=True):
        assert point in self
        bs = point.bs
        rv = (
            LatticeSite((bs[0], bs[1] + 1), 0),
            LatticeSite((bs[0] + 1, bs[1] + 1), 0),
            LatticeSite((bs[0] + 1, bs[1]), 0),
        )
        if double_count:
            rv += (
                LatticeSite((bs[0], bs[1] - 1), 0),
                LatticeSite((bs[0] - 1, bs[1] - 1), 0),
                LatticeSite((bs[0] - 1, bs[1]), 0),
            )
        return rv

    def next_nearest_neighbors(self, point, double_count=True):
        assert point in self
        bs = point.bs
        rv = (
            LatticeSite((bs[0] + 2, bs[1] + 1), 0),
            LatticeSite((bs[0] + 1, bs[1] + 2), 0),
            LatticeSite((bs[0] - 1, bs[1] + 1), 0),
        )
        if double_count:
            rv += (
                LatticeSite((bs[0] - 2, bs[1] - 1), 0),
                LatticeSite((bs[0] - 1, bs[1] - 2), 0),
                LatticeSite((bs[0] + 1, bs[1] - 1), 0),
            )
        return rv
