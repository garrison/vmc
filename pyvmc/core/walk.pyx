import abc

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.utils.immutable import Immutable

class WalkPlan(Immutable):
    __slots__ = ("wavefunction",)

    def init_validate(self, wavefunction, *args, **kwargs):
        assert isinstance(wavefunction, Wavefunction)
        return super(WalkPlan, self).init_validate(wavefunction, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

class StandardWalkPlan(WalkPlan):
    __slots__ = ("wavefunction",)

    def to_json(self):
        return {"walk-type": "standard"}
