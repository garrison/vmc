"""Classes representing a spatial subset of a lattice.

Often used e.g. in Renyi entropy calculations.
"""

import abc
import numbers
import collections

from cython.operator cimport dereference as deref

import numpy

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

    def __hash__(self):
        return hash(tuple(self)) | hash(self.lattice)

    def __richcmp__(self, other, int op):
        if op == 2:  # ==
            return (isinstance(other, Subsystem) and
                    tuple(self) == tuple(other) and
                    self.lattice == other.lattice)
        elif op == 3:  # !=
            return (not isinstance(other, Subsystem) or
                    tuple(self) != tuple(other) or
                    self.lattice != other.lattice)
        # we don't implement <, <=, >, >=
        raise NotImplementedError

# NOTE: each subclass of Subsystem should implement all the Hashable and
# Sequence methods.
collections.Hashable.register(Subsystem)
collections.Sequence.register(Subsystem)

cdef class SimpleSubsystem(Subsystem):
    """Subsystem consisting of a hyper-rectangle bordering the origin"""

    def __init__(self, dimensions, lattice):
        super(SimpleSubsystem, self).__init__(lattice)
        assert isinstance(dimensions, collections.Sequence)
        if len(dimensions) > MAX_DIMENSION:
            raise ValueError("provided subsystem has greater than {} dimensions".format(MAX_DIMENSION))
        assert all(isinstance(d, numbers.Integral) and d > 0 for d in dimensions)
        if any(d1 > d2 for d1, d2 in zip(dimensions, lattice.dimensions)):
            raise Exception("subsystem cannot be larger than the system")
        cdef CppSimpleSubsystemDimensionVector v
        for i, x in enumerate(dimensions):
            v.push_back(x)
        self.sharedptr.reset(new CppSimpleSubsystem(v))

    property dimensions:
        def __get__(self):
            cdef const CppSimpleSubsystemDimensionVector *v = &(<CppSimpleSubsystem*>self.sharedptr.get()).subsystem_length
            cdef unsigned int i
            return tuple([deref(v)[i] for i in xrange(v.size())])

    def __len__(self):
        return numpy.product(self.dimensions) * self.lattice.basis_indices

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

cdef class CustomSubsystem(Subsystem):
    """Subsystem built from a function that tests site membership

    >>> subsystem = CustomSubsystem(lambda site: sum(site.bs) % 2 == 0, Lattice([4, 4]))
    """

    cdef tuple _sites
    cdef tuple _site_indices

    def __init__(self, contains_lambda, lattice):
        super(CustomSubsystem, self).__init__(lattice)

        cdef unsigned int i
        cdef bint b
        cdef dynamic_bitset site_status_array

        site_status_array.resize(len(lattice))
        _sites = []
        _site_indices = []

        for i, site in enumerate(lattice):
            b = bool(contains_lambda(site))
            site_status_array[i] = b
            if b:
                _sites.append(lattice[i])
                _site_indices.append(i)

        self.sharedptr.reset(new CppCustomSubsystem(site_status_array))
        self._sites = tuple(_sites)
        self._site_indices = tuple(_site_indices)

    def __len__(self):
        return len(self._sites)

    def __iter__(self):
        return iter(self._sites)

    def __getitem__(self, index):
        return self._sites[index]

    def index(self, site):
        return self._sites.index(site)

    def __contains__(self, site):
        return site in self._sites

    def count(self, site):
        return 1 if site in self else 0

    def __repr__(self):
        return "{}(lambda s: {lattice}.index(s) in {site_indices}, {lattice})".format(self.__class__.__name__,
                                   site_indices=repr(self._site_indices),
                                   lattice=repr(self.lattice))

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
            ("_site_indices", self._site_indices),
        ])

    @staticmethod
    def _from_json(json_repr, lattice):
        assert json_repr["type"] == "CustomSubsystem"
        return CustomSubsystem(lambda s: lattice.index(s) in json_repr["_site_indices"], lattice)

__subsystem_registry["CustomSubsystem"] = CustomSubsystem
