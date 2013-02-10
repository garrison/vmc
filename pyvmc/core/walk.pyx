from pyvmc.includes.boost.shared_ptr cimport shared_ptr

import abc
import collections

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.wavefunction cimport CppWavefunctionAmplitude, create_wfa
from pyvmc.core.rng cimport RandomNumberGenerator
from pyvmc.utils.immutable import Immutable

class WalkPlan(Immutable):
    __slots__ = ("wavefunction",)

    def init_validate(self, wavefunction, *args, **kwargs):
        assert isinstance(wavefunction, Wavefunction)
        return super(WalkPlan, self).init_validate(wavefunction, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
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

    def create_walk(self, RandomNumberGenerator rng not None):
        cdef shared_ptr[CppWavefunctionAmplitude] wfa = create_wfa(self.wavefunction, rng)
        cdef Walk walk = Walk()
        walk.autoptr.reset(new CppStandardWalk(wfa))
        return walk
