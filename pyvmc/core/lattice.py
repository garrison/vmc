"""
"""

import abc
import numbers
import collections

from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.utils import product

class BravaisSite(tuple):
    pass

_LatticeSite = collections.namedtuple('LatticeSite', ('bs', 'bi'))

class LatticeSite(_LatticeSite):
    """represents a site on a lattice

    bs == bravais site
    bi == basis index
    """

    __slots__ = ()

    def __new__(cls, bs, bi=0):
        assert isinstance(bs, collections.Sequence)
        bs = BravaisSite(bs)
        assert all(isinstance(x, numbers.Integral) for x in bs)
        assert isinstance(bi, numbers.Integral) and bi >= 0
        return _LatticeSite.__new__(cls, bs, bi)

class Lattice(collections.Sequence, collections.Hashable):
    __slots__ = ('dimensions', 'basis_indices')

    def __init__(self, dimensions, basis_indices=1):
        assert isinstance(dimensions, collections.Sequence)
        assert len(dimensions) != 0
        dimensions = tuple(dimensions)
        assert all(isinstance(x, numbers.Integral) and x > 0 for x in dimensions)
        object.__setattr__(self, "dimensions", dimensions)

        assert isinstance(basis_indices, numbers.Integral) and basis_indices > 0
        object.__setattr__(self, "basis_indices", basis_indices)

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

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and
                self.dimensions == other.dimensions and
                self.basis_indices == other.basis_indices)

    def __ne__(self, other):
        return (self is not other) and not self.__eq__(other)

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

    def __setattr__(self, name, value):
        raise TypeError

    def __delattr__(self, name):
        raise TypeError

    def __hash__(self):
        return hash(self.dimensions) | hash(self.basis_indices)

    def __repr__(self):
        return "%s(%s, %s)" % (
            self.__class__.__name__,
            repr(self.dimensions),
            repr(self.basis_indices)
        )

class LatticeRealization(Lattice):
    __metaclass__ = abc.ABCMeta

    __slots__ = ('dimensions', 'basis_indices')

    def abstract_lattice(self):
        return Lattice(self.dimensions, self.basis_indices)

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

    def __init__(self, dimensions):
        super(HypercubicLattice, self).__init__(dimensions, 1)

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

    def __init__(self, dimensions):
        assert len(dimensions) == 2
        super(HexagonalLattice, self).__init__(dimensions, 1)

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
