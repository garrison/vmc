from pyvmc.core.measurement import WalkPlan, MeasurementPlan
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.subsystem import Subsystem

class RenyiModWalkPlan(WalkPlan):
    # NOTE: this is deprecated in favor of renyi-mod/possible

    __slots__ = ("wavefunction",)

    def to_json(self):
        return {"walk-type": "renyi-mod"}

class RenyiModMeasurementPlan(MeasurementPlan):
    # NOTE: this is deprecated in favor of renyi-mod/possible

    __slots__ = ("walk", "subsystem")

    def __init__(self, wavefunction, subsystem):
        walk = RenyiModWalkPlan(wavefunction)
        super(RenyiModMeasurementPlan, self).__init__(walk, subsystem)

class RenyiModPossibleWalkPlan(WalkPlan):
    __slots__ = ("wavefunction", "subsystem")

    def init_validate(self, wavefunction, subsystem):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(subsystem, Subsystem)
        assert wavefunction.lattice == subsystem.lattice
        return wavefunction, subsystem

    def to_json(self):
        return {
            "walk-type": "renyi-mod/possible",
            "subsystem": self.subsystem.to_json(),
        }

class RenyiModPossibleMeasurementPlan(MeasurementPlan):
    __slots__ = ("walk",)

    def __init__(self, wavefunction, subsystem):
        walk = RenyiModPossibleWalkPlan(wavefunction, subsystem)
        super(RenyiModPossibleMeasurementPlan, self).__init__(walk)

    def to_json(self):
        return {"type": "renyi-mod/possible"}

class RenyiSignWalkPlan(WalkPlan):
    __slots__ = ("wavefunction", "subsystem")

    def init_validate(self, wavefunction, subsystem):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(subsystem, Subsystem)
        assert wavefunction.lattice == subsystem.lattice
        return wavefunction, subsystem

    def to_json(self):
        return {
            "walk-type": "renyi-sign",
            "subsystem": self.subsystem.to_json(),
        }

class RenyiSignMeasurementPlan(MeasurementPlan):
    __slots__ = ("walk",)

    def __init__(self, wavefunction, subsystem):
        walk = RenyiSignWalkPlan(wavefunction, subsystem)
        super(RenyiSignMeasurementPlan, self).__init__(walk)

    def to_json(self):
        return {"type": "renyi-sign"}
