"""
"""

import abc
import numbers
import collections

from pyvmc.utils import product
from pyvmc.core.lattice import LatticeSite, Lattice
from pyvmc.utils.immutable import Immutable

class Subsystem(Immutable, collections.Sequence):
    """Abstract base class representing a spatial subset of a lattice"""

    __slots__ = ('lattice',)

    def init_validate(self, lattice):
        assert isinstance(lattice, Lattice)
        return (lattice,)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

    def __len__(self):
        count = 0
        for site in self.lattice:
            if site in self:
                count += 1
        return count

class SimpleSubsystem(Subsystem):
    """Subsystem consisting of a hyper-rectangle bordering the origin"""

    __slots__ = ('dimensions', 'lattice')

    def init_validate(self, dimensions, lattice):
        (lattice,) = super(SimpleSubsystem, self).init_validate(lattice)
        assert isinstance(dimensions, collections.Sequence)
        assert all(isinstance(d, numbers.Integral) and d > 0 for d in dimensions)
        if any(d1 > d2 for d1, d2 in zip(dimensions, lattice.dimensions)):
            raise Exception("subsystem cannot be larger than the system")
        return tuple(dimensions), lattice

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
        index_offset = 1
        for i, x in enumerate(site.bs):
            index += x * index_offset
            index_offset *= dimensions[i]
        index += site.bi * index_offset
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
