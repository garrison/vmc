"""All the Lattice infrastructure is in this file.

This wraps Lattice.hpp and provides abstractions on top of it, e.g. various LatticeRealization's.
"""

import six

import abc
import numbers
import collections

import numpy
from cython.operator cimport dereference as deref

from pyvmc.core.boundary_conditions import valid_boundary_conditions, periodic_bc, antiperiodic_bc, open_bc
from pyvmc.utils import product
from pyvmc.constants import two_pi, sqrt_three_over_two, pi

cdef class LatticeSite(object):
    """represents a site on a lattice

    bs == bravais site
    bi == basis index

    remains constant after initialization
    """

    def __init__(self, bs, bi=0):
        assert isinstance(bs, collections.Sequence)
        assert len(bs) > 0
        assert all(isinstance(x, numbers.Integral) for x in bs)
        assert isinstance(bi, numbers.Integral) and bi >= 0
        if len(bs) > MAX_DIMENSION:
            raise ValueError("provided site has greater than {} dimensions".format(MAX_DIMENSION))
        self.cpp.set_n_dimensions(len(bs))
        cdef unsigned int i
        cdef int x
        for i, x in enumerate(bs):
            self.cpp.set_bs_coordinate(i, x)
        self.cpp.basis_index = bi

    property bs:
        def __get__(self):
            cdef unsigned int i
            return tuple([self.cpp[i] for i in range(self.cpp.n_dimensions())])

    property bi:
        def __get__(self):
            return self.cpp.basis_index

    def displaced(self, displacement):
        """Returns the site displaced by a Bravais lattice vector (given by a tuple of integers)"""
        assert isinstance(displacement, collections.Sequence)
        assert all(isinstance(x, numbers.Integral) for x in displacement)
        if len(displacement) != self.cpp.n_dimensions():
            raise ValueError("Displacement vector has wrong number of dimensions")
        return LatticeSite([a + b for a, b in zip(self.bs, displacement)],
                           self.bi)

    def negative_displaced(self, displacement):
        """Returns the site displaced in the negative direction by a Bravais lattice vector (given by a tuple of integers)"""
        assert isinstance(displacement, collections.Sequence)
        assert all(isinstance(x, numbers.Integral) for x in displacement)
        if len(displacement) != self.cpp.n_dimensions():
            raise ValueError("Displacement vector has wrong number of dimensions")
        return LatticeSite([a - b for a, b in zip(self.bs, displacement)],
                           self.bi)

    def __add__(self, displacement):
        """Adds a displacement to the given LatticeSite and returns the result.

        The basis index will remain the same.

        The LatticeSite must be given on the left side of the expression.

        >>> LatticeSite([1, 2]) + (2, 3)
        LatticeSite((3, 5), 0)
        """
        if not isinstance(self, LatticeSite):
            # see http://docs.cython.org/src/userguide/special_methods.html#arithmetic-methods
            return NotImplemented
        return self.displaced(displacement)

    def __sub__(self, displacement):
        """Subtracts a displacement from the given LatticeSite and returns the result.

        The basis index will remain the same.

        The LatticeSite must be given on the left side of the expression.

        >>> LatticeSite([2, 4], 1) - [1, -3]
        LatticeSite((1, 7), 1)
        """
        if not isinstance(self, LatticeSite):
            # see http://docs.cython.org/src/userguide/special_methods.html#arithmetic-methods
            return NotImplemented
        return self.negative_displaced(displacement)

    def __hash__(self):
        return hash(self.bs) | hash(self.bi)

    def __richcmp__(LatticeSite self not None, other, int op):
        if self.__class__ != other.__class__:
            if op == 2:  #  ==
                return False
            elif op == 3: #  !=
                return True
            else:
                raise NotImplementedError

        cdef LatticeSite other_ = other

        if op == 2:  # ==
            return self.cpp == other_.cpp
        elif op == 3:  # !=
            return self.cpp != other_.cpp
        if op == 0:  # <
            return self.cpp < other_.cpp
        elif op == 1:  # <=
            return not other_.cpp < self.cpp
        elif op == 4:  # >
            return other_.cpp < self.cpp
        elif op == 5:  # >=
            return not self.cpp < other_.cpp

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.bs),
                                   self.bi)

    def to_json(self):
        return [self.bs, self.bi]

cdef LatticeSite_from_cpp(CppLatticeSite cpp_lattice_site):
    cdef LatticeSite lattice_site = LatticeSite.__new__(LatticeSite)
    lattice_site.cpp = cpp_lattice_site
    return lattice_site

collections.Hashable.register(LatticeSite)

cdef __phase_wrap(wraps, phase):
    if wraps == 0:
        return 1
    elif phase == 0:
        # consider this case explicitly so we never raise 0 to a negative power
        return 0
    else:
        return phase ** wraps

cdef class Lattice(object):
    def __init__(self, dimensions, basis_indices=1):
        assert isinstance(dimensions, collections.Sequence)
        assert len(dimensions) != 0
        assert all(isinstance(x, numbers.Integral) and x > 0 for x in dimensions)
        assert isinstance(basis_indices, numbers.Integral) and basis_indices > 0
        cdef DimensionVector v
        for i, x in enumerate(dimensions):
            v.push_back(x)
        self.sharedptr.reset(new CppLattice(v, basis_indices))

    property dimensions:
        def __get__(self):
            cdef unsigned int i
            return tuple([self.sharedptr.get().dimensions[i] for i in range(self.sharedptr.get().n_dimensions())])

    property basis_indices:
        def __get__(self):
            return self.sharedptr.get().basis_indices

    def to_json(self):
        return collections.OrderedDict([
            ("dimensions", list(self.dimensions)),
            ("basis_indices", self.basis_indices),
        ])

    @classmethod
    def from_json(cls, d):
        assert cls == Lattice  # for now
        return Lattice(d["dimensions"], d["basis_indices"])

    def abstract_lattice(self):
        """This will always return a Lattice object, never a LatticeRealization"""
        return self

    def enforce_boundary(self, site, boundary_conditions):
        """Enforce the boundary of a site which may be outside the lattice

        Any quantum amplitudes should be multiplied by `phase`
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
        if boundary_conditions is False:
            return new_site
        assert valid_boundary_conditions(boundary_conditions, len(lattice_dimensions))
        phase = product(__phase_wrap(x // length, bc.phase) for x, length, bc
                        in zip(bravais_site, lattice_dimensions, boundary_conditions))
        return new_site, phase

    def __len__(self):
        return self.sharedptr.get().total_sites()

    def total_bravais_sites(self):
        return self.sharedptr.get().total_bravais_sites()

    def __iter__(self):
        for bi in xrange(self.basis_indices):
            for bs in self.iterate_bravais_sites():
                yield LatticeSite(bs, bi)

    def iterate_bravais_sites(self):
        dimensions = self.dimensions
        d = len(dimensions)
        s = [0] * d
        while True:
            yield tuple(s)
            for i in xrange(d):
                s[i] += 1
                if s[i] < dimensions[i]:
                    break
                s[i] = 0
            else:
                raise StopIteration

    def __getitem__(self, unsigned int index):
        if index >= len(self):
            raise ValueError
        return LatticeSite_from_cpp(self.sharedptr.get().site_from_index(index))

    def index(self, LatticeSite site not None):
        if site not in self:
            raise ValueError
        return <int>self.sharedptr.get().site_to_index(site.cpp)

    def __contains__(self, LatticeSite site not None):
        return bool(site.cpp.n_dimensions() == self.sharedptr.get().n_dimensions() and
                    self.sharedptr.get().site_is_valid(site.cpp))

    def count(self, x):
        return 1 if x in self else 0

    def __hash__(self):
        return hash(self.dimensions) | hash(self.basis_indices)

    def __richcmp__(self, other, int op):
        if op == 2:  # ==
            return (self.__class__ == other.__class__ and
                    self.dimensions == other.dimensions and
                    self.basis_indices == other.basis_indices)
        elif op == 3:  # !=
            return (self.__class__ != other.__class__ or
                    self.dimensions != other.dimensions or
                    self.basis_indices != other.basis_indices)
        # we don't implement <, <=, >, >=
        raise NotImplementedError

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.dimensions),
                                   self.basis_indices)

collections.Hashable.register(Lattice)
collections.Sequence.register(Lattice)

class LatticeRealization(six.with_metaclass(abc.ABCMeta, Lattice)):
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

    def spatial_coordinates(self, point):
        """
        Returns the real-space coordinates of the given LatticeSite

        NOTE: this method doesn't even require that the resulting point be on
        the lattice
        """
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
    a = 1.0   # lattice spacing

    def __init__(self, dimensions):
        super(HypercubicLattice, self).__init__(dimensions, 1)

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
    a = 1.0   # lattice spacing

    def __init__(self, dimensions):
        assert len(dimensions) == 2
        super(HexagonalLattice, self).__init__(dimensions, 1)

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
            LatticeSite((bs[0] + 1, bs[1]), 0),
            LatticeSite((bs[0] + 1, bs[1] + 1), 0),
            LatticeSite((bs[0], bs[1] + 1), 0),
        )
        if double_count:
            rv += (
                LatticeSite((bs[0] - 1, bs[1]), 0),
                LatticeSite((bs[0] - 1, bs[1] - 1), 0),
                LatticeSite((bs[0], bs[1] - 1), 0),
            )
        return rv

    def second_nearest_neighbors(self, point, double_count=True):
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

    next_nearest_neighbors = second_nearest_neighbors

    def third_nearest_neighbors(self, point, double_count=True):
        assert point in self
        bs = point.bs
        rv = (
            LatticeSite((bs[0] + 2, bs[1]), 0),
            LatticeSite((bs[0] + 2, bs[1] + 2), 0),
            LatticeSite((bs[0], bs[1] + 2), 0),
        )
        if double_count:
            rv += (
                LatticeSite((bs[0] - 2, bs[1]), 0),
                LatticeSite((bs[0] - 2, bs[1] - 2), 0),
                LatticeSite((bs[0], bs[1] - 2), 0),
            )
        return rv

    def basic_plaquettes(self):
        return (
            (LatticeSite([0, 0]), LatticeSite([1, 0]), LatticeSite([1, 1]), LatticeSite([0, 1])),
            (LatticeSite([0, 0]), LatticeSite([1, 0]), LatticeSite([2, 1]), LatticeSite([1, 1])),
            (LatticeSite([0, 0]), LatticeSite([1, 1]), LatticeSite([1, 2]), LatticeSite([0, 1])),
        )
