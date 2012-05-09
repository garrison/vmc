"""
"""

import abc
import numbers
import collections

from pyvmc.utils import product
from pyvmc.core.lattice import LatticeSite, Lattice

class Subsystem(collections.Sequence, collections.Hashable):
    """Abstract base class representing a spatial subset of a lattice"""

    __metaclass__ = abc.ABCMeta
    __slots__ = ('lattice',)

    def __init__(self, lattice):
        assert isinstance(lattice, Lattice)
        object.__setattr__(self, 'lattice', lattice)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

    def __len__(self):
        count = 0
        for site in self.lattice:
            if site in self:
                count += 1
        return count

    def __setattr__(self, name, value):
        raise TypeError

    def __delattr__(self, name):
        raise TypeError

class SimpleSubsystem(Subsystem):
    """Subsystem consisting of a hyper-rectangle bordering the origin"""

    __slots__ = ('dimensions', 'lattice')

    def __init__(self, dimensions, lattice):
        super(SimpleSubsystem, self).__init__(lattice)
        object.__setattr__(self, 'dimensions', tuple(dimensions))
        assert all(isinstance(d, numbers.Integral) and d > 0 for d in self.dimensions)
        if any(d1 > d2 for d1, d2 in zip(dimensions, lattice.dimensions)):
            raise Exception("subsystem cannot be larger than the system")

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and
                self.dimensions == other.dimensions and
                self.lattice == other.lattice)

    def __ne__(self, other):
        return (self is not other) and not self.__eq__(other)

    def __len__(self):
        return product(self.dimensions) * self.lattice.basis_indices

    def __iter__(self):
        dimensions = self.dimensions
        d = len(dimensions)
        s_bound = dimensions + (self.lattice.basis_indices,)
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
        return bool(site in self.lattice and
                    all(x < y for x, y in zip(site.bs, self.dimensions)))

    def to_json(self):
        return {
            "type": "simple",
            "dimensions": self.dimensions,
        }

    def __hash__(self):
        return hash(self.dimensions) | hash(self.lattice)

    def __repr__(self):
        return "%s(%s, %s)" % (
            self.__class__.__name__,
            repr(self.dimensions),
            repr(self.lattice)
        )
