import six

import abc
import numbers
import collections
from itertools import chain

from pyvmc.utils.immutable import Immutable
from pyvmc.core.lattice import Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.wavefunction import Wavefunction

class SiteHop(Immutable):
    __slots__ = ("source", "destination", "species")

    def init_validate(self, source, destination, species):
        assert isinstance(source, LatticeSite)
        assert isinstance(destination, LatticeSite)
        assert isinstance(species, numbers.Integral) and species >= 0
        return (source, destination, species)

    def is_valid_for(self, wavefunction):
        assert isinstance(wavefunction, Wavefunction)
        return (self.source in wavefunction.lattice and
                self.destination in wavefunction.lattice and
                self.species < wavefunction.N_species)

    def to_json(self):
        return collections.OrderedDict([
            ("source", self.source.to_json()),
            ("destination", self.destination.to_json()),
            ("species", self.species),
        ])

class Operator(six.with_metaclass(abc.ABCMeta)):
    @abc.abstractmethod
    def get_basic_operators(self):
        raise NotImplementedError

    @abc.abstractmethod
    def evaluate(self, context):
        """context is a nested dict where a BasicOperator is the key and the value is the expectation value of that operator"""
        raise NotImplementedError

class BasicOperator(Immutable):
    """A BasicOperator represents anything that can be represented a SiteHop's

    If a sum is to be performed over all sites, set boundary_conditions.
    Otherwise, if no sum is to be performed set boundary_conditions to None.
    """

    __slots__ = ("hops", "boundary_conditions")

    def init_validate(self, hops, boundary_conditions=None):
        assert isinstance(hops, collections.Sequence)
        hops = tuple(sorted(hops, key=lambda hop: (hop.species, hop.source, hop.destination)))
        assert all([isinstance(hop, SiteHop) for hop in hops])
        if boundary_conditions is not None:
            boundary_conditions = tuple(boundary_conditions)
            assert valid_boundary_conditions(boundary_conditions, len(boundary_conditions))
        return (hops, boundary_conditions)

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
            ("hops", [hop.to_json() for hop in self.hops]),
            ("boundary_conditions", self.boundary_conditions),
        ])

    # these next two methods exist only so we can pass around Operator's
    # without caring whether they are BasicOperator's or CompositeOperator's
    # (see the Operator abstact base class, which declares these methods).

    def get_basic_operators(self):
        return {self}

    def evaluate(self, context):
        def _evaluate():
            return context[self]
        return _evaluate

Operator.register(BasicOperator)

class CompositeOperator(Immutable):
    """
    A CompositeOperator represents anything that is built from BasicOperators
    (e.g. a linear combination of them)
    """

    __slots__ = ("operators",)

    parameters = ()

    def get_basic_operators(self):
        return set(chain.from_iterable([o.get_basic_operators() for o in self.operators]))

Operator.register(CompositeOperator)
