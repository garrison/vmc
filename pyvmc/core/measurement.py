import abc
import collections

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.subsystem import Subsystem
from pyvmc.utils.immutable import Immutable

class WalkPlan(Immutable):
    __slots__ = ("wavefunction",)

    def init_validate(self, wavefunction, *args, **kwargs):
        assert isinstance(wavefunction, Wavefunction)
        return super(WalkPlan, self).init_validate(wavefunction, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

class MeasurementPlan(Immutable):
    __slots__ = ("walk",)

    def init_validate(self, walk, *args, **kwargs):
        assert isinstance(walk, WalkPlan)
        return super(MeasurementPlan, self).init_validate(walk, *args, **kwargs)

    @abc.abstractmethod
    def to_json(self):
        raise NotImplementedError

class StandardWalkPlan(WalkPlan):
    __slots__ = ("wavefunction",)

    def to_json(self):
        return {"walk-type": "standard"}

class DensityDensityMeasurementPlan(MeasurementPlan):
    __slots__ = ()

class GreenMeasurementPlan(MeasurementPlan):
    __slots__ = ()

class SubsystemOccupationProbabilityMeasurementPlan(MeasurementPlan):
    __slots__ = ("walk", "subsystem")

    def __init__(self, wavefunction, subsystem):
        walk = StandardWalkPlan(wavefunction)
        super(SubsystemOccupationProbabilityMeasurementPlan, self).__init__(walk, subsystem)

    def to_json(self):
        return {
            "type": "subsystem-occupation-number-probability",
            "subsystem": self.subsystem.to_json(),
            "steps-per-measurement": 100
        }
