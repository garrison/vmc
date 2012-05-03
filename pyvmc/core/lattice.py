"""
"""

import abc
import numbers
import collections

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
