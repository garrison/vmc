import six

from pyvmc.includes.boost.shared_ptr cimport shared_ptr

import abc
import collections

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.wavefunction cimport CppWavefunctionAmplitude, create_wfa
from pyvmc.core.rng cimport RandomNumberGenerator
from pyvmc.utils.immutable import Immutable, ImmutableMetaclass

__walk_plan_registry = {}

class WalkPlanMetaclass(ImmutableMetaclass):
    def __init__(cls, name, bases, dct):
        assert name not in __walk_plan_registry
        __walk_plan_registry[name] = cls
        super(WalkPlanMetaclass, cls).__init__(name, bases, dct)

class WalkPlan(six.with_metaclass(WalkPlanMetaclass, Immutable)):
    __slots__ = ("wavefunction",)

    def init_validate(self, wavefunction, *args, **kwargs):
        assert isinstance(wavefunction, Wavefunction)
        return super(WalkPlan, self).init_validate(wavefunction, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

    @classmethod
    def from_json(cls, json_repr, wavefunction):
        cls_ = __walk_plan_registry[json_repr["type"]]
        assert issubclass(cls_, cls)
        return cls_._from_json(json_repr, wavefunction)

    @staticmethod
    def _from_json(json_repr, wavefunction):
        raise NotImplementedError

    @abc.abstractmethod
    def create_walk(self, RandomNumberGenerator rng not None):
        raise NotImplementedError

class StandardWalkPlan(WalkPlan):
    __slots__ = ("wavefunction",)

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
        ])

    @staticmethod
    def _from_json(json_repr, wavefunction):
        assert json_repr["type"] == "StandardWalkPlan"
        return StandardWalkPlan(wavefunction)

    def create_walk(self, RandomNumberGenerator rng not None):
        cdef shared_ptr[CppWavefunctionAmplitude] wfa = create_wfa(self.wavefunction, rng)
        cdef Walk walk = Walk()
        walk.autoptr.reset(new CppStandardWalk(wfa))
        return walk
