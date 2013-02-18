"""
"""

import abc
import numbers
import collections

from cython.operator cimport dereference as deref

from pyvmc.utils import product
from pyvmc.core.lattice import LatticeSite, Lattice
from pyvmc.utils.immutable import Immutable

__subsystem_registry = {}

# commented out because metaclasses don't seem to work with cython.  for now we
# explicitly register the subclasses.
#
#class SubsystemMetaclass(type):
#    def __init__(cls, name, bases, dct):
#        assert name not in __subsystem_registry
#        __subsystem_registry[name] = cls
#        super(SubsystemMetaclass, cls).__init__(name, bases, dct)

cdef class Subsystem(object):
    """Abstract base class representing a spatial subset of a lattice"""

    #__metaclass__ = SubsystemMetaclass

    #__metaclass__ = abc.ABCMeta

    def __init__(self, lattice):
        assert isinstance(lattice, Lattice)
        self.lattice_ = lattice

    property lattice:
        def __get__(self):
            return self.lattice_

    # The following line is commented out because it does not work with Cython.
    #@abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

    @classmethod
    def from_json(cls, json_repr, lattice):
        cls_ = __subsystem_registry[json_repr["type"]]
        assert issubclass(cls_, cls)
        return cls_._from_json(json_repr, lattice)

    def __len__(self):
        cdef int count = 0
        for site in self.lattice:
            if site in self:
                count += 1
        return count

# NOTE: each subclass of Subsystem should implement all the Hashable and
# Sequence methods.
collections.Hashable.register(Subsystem)
collections.Sequence.register(Subsystem)

cdef class SimpleSubsystem(Subsystem):
    """Subsystem consisting of a hyper-rectangle bordering the origin"""

    def __init__(self, dimensions, lattice):
        super(SimpleSubsystem, self).__init__(lattice)
        assert isinstance(dimensions, collections.Sequence)
        assert all(isinstance(d, numbers.Integral) and d > 0 for d in dimensions)
        if any(d1 > d2 for d1, d2 in zip(dimensions, lattice.dimensions)):
            raise Exception("subsystem cannot be larger than the system")
        cdef UDimensionVector v
        for i, x in enumerate(dimensions):
            v.push_back(x)
        self.sharedptr.reset(new CppSimpleSubsystem(v))

    property dimensions:
        def __get__(self):
            cdef const UDimensionVector *v = &(<CppSimpleSubsystem*>self.sharedptr.get()).subsystem_length
            cdef unsigned int i
            return tuple([deref(v)[i] for i in xrange(v.size())])

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

    def count(self, site):
        return 1 if site in self else 0

    def __hash__(self):
        return hash(self.dimensions) | hash(self.lattice)

    def __richcmp__(self, other, int op):
        if op == 2:  # ==
            return (self.__class__ == other.__class__ and
                    self.dimensions == other.dimensions and
                    self.lattice == other.lattice)
        elif op == 3:  # !=
            return (self.__class__ != other.__class__ or
                    self.dimensions != other.dimensions or
                    self.lattice != other.lattice)
        # we don't implement <, <=, >, >=
        raise NotImplementedError

    def __repr__(self):
        return "{}({}, {})".format(self.__class__.__name__,
                                   repr(self.dimensions),
                                   repr(self.lattice))

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
            ("dimensions", self.dimensions),
        ])

    @staticmethod
    def _from_json(json_repr, lattice):
        assert json_repr["type"] == "SimpleSubsystem"
        return SimpleSubsystem(json_repr["dimensions"], lattice)

__subsystem_registry["SimpleSubsystem"] = SimpleSubsystem
