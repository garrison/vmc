from pyvmc.core.measurement import WalkPlan, MeasurementPlan
from pyvmc.core.measurement cimport BaseMeasurement
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.subsystem cimport Subsystem
from pyvmc.core cimport complex_t

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

    def to_measurement(self):
        return RenyiModPossibleMeasurement()

cdef class RenyiModPossibleMeasurement(BaseMeasurement):
    def __init__(self):
        self.sharedptr.reset(new CppRenyiModPossibleMeasurement())

    def get_result(self):
        return (<CppRenyiModPossibleMeasurement*>self.sharedptr.get()).get_estimate().get_result()

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

    def to_measurement(self):
        return RenyiSignMeasurement()

cdef class RenyiSignMeasurement(BaseMeasurement):
    def __init__(self):
        self.sharedptr.reset(new CppRenyiSignMeasurement())

    def get_result(self):
        cdef complex_t c = (<CppRenyiSignMeasurement*>self.sharedptr.get()).get_estimate().get_result()
        return complex(c.real(), c.imag())
