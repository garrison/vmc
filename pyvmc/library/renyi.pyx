from cython.operator cimport dereference as deref

from collections import OrderedDict

from pyvmc.core.walk import WalkPlan
from pyvmc.core.walk cimport Walk
from pyvmc.core.measurement import BasicMeasurementPlan, CompositeMeasurementPlan
from pyvmc.core.measurement cimport BaseMeasurement
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.wavefunction cimport CppWavefunctionAmplitude, WavefunctionWrapper, std_move_wfa
from pyvmc.core.estimate cimport Estimate_from_CppRealBlockedEstimate, Estimate_from_CppComplexBlockedEstimate
from pyvmc.core.subsystem cimport Subsystem
from pyvmc.core.rng cimport RandomNumberGenerator
from pyvmc.core cimport complex_t
from pyvmc.measurements import SubsystemOccupationProbabilityMeasurementPlan

class RenyiModPossibleWalkPlan(WalkPlan):
    __slots__ = ("wavefunction", "subsystem")

    def init_validate(self, wavefunction, subsystem):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(subsystem, Subsystem)
        assert wavefunction.lattice == subsystem.lattice
        return wavefunction, subsystem

    def to_json(self):
        return OrderedDict([
            ("type", self.__class__.__name__),
            ("subsystem", self.subsystem.to_json()),
        ])

    @staticmethod
    def _from_json(json_repr, wavefunction):
        assert json_repr["type"] == "RenyiModPossibleWalkPlan"
        return RenyiModPossibleWalkPlan(wavefunction, Subsystem.from_json(json_repr["subsystem"], wavefunction.lattice))

    def create_walk(self, RandomNumberGenerator rng not None):
        assert rng.is_good()
        cdef Subsystem subsystem = self.subsystem
        cdef WavefunctionWrapper ww = self.wavefunction.to_wavefunction()
        cdef Walk walk = Walk()
        # We need two copies of the system, each of which has the same number
        # of particles in the subsystem.  So for now we just initialize both
        # copies with the same exact positions.
        walk.uniqueptr.reset(new CppRenyiModPossibleWalk(std_move_wfa(ww.sharedptr.get().create_nonzero_wavefunctionamplitude(ww.sharedptr, deref(rng.uniqueptr.get()))), subsystem.sharedptr))
        return walk

class RenyiModPossibleMeasurementPlan(BasicMeasurementPlan):
    __slots__ = ("walk_plan",)

    def __init__(self, wavefunction, subsystem):
        walk_plan = RenyiModPossibleWalkPlan(wavefunction, subsystem)
        super(RenyiModPossibleMeasurementPlan, self).__init__(walk_plan)

    def to_json(self):
        return OrderedDict([
            ("type", self.__class__.__name__),
            ("subsystem", self.walk_plan.subsystem.to_json()),
        ])

    @staticmethod
    def _from_json(json_repr, wavefunction):
        assert json_repr["type"] == "RenyiModPossibleMeasurementPlan"
        return RenyiModPossibleMeasurementPlan(wavefunction, Subsystem.from_json(json_repr["subsystem"], wavefunction.lattice))

    def to_measurement(self):
        return RenyiModPossibleMeasurement()

cdef class RenyiModPossibleMeasurement(BaseMeasurement):
    def __init__(self):
        self.sharedptr.reset(new CppRenyiModPossibleMeasurement())

    def get_estimate(self, key=None):
        if key is not None:
            raise KeyError
        return Estimate_from_CppRealBlockedEstimate((<CppRenyiModPossibleMeasurement*>self.sharedptr.get()).get_estimate())

    def get_estimates(self):
        return {None: self.get_estimate()}

class RenyiSignWalkPlan(WalkPlan):
    __slots__ = ("wavefunction", "subsystem")

    def init_validate(self, wavefunction, subsystem):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(subsystem, Subsystem)
        assert wavefunction.lattice == subsystem.lattice
        return wavefunction, subsystem

    def to_json(self):
        return OrderedDict([
            ("type", self.__class__.__name__),
            ("subsystem", self.subsystem.to_json()),
        ])

    @staticmethod
    def _from_json(json_repr, wavefunction):
        assert json_repr["type"] == "RenyiSignWalkPlan"
        return RenyiSignWalkPlan(wavefunction, Subsystem.from_json(json_repr["subsystem"], wavefunction.lattice))

    def create_walk(self, RandomNumberGenerator rng not None):
        assert rng.is_good()
        cdef Subsystem subsystem = self.subsystem
        cdef WavefunctionWrapper ww = self.wavefunction.to_wavefunction()
        cdef Walk walk = Walk()
        # We need two copies of the system, each of which has the same number
        # of particles in the subsystem.  So for now we just initialize both
        # copies with the same exact positions.
        walk.uniqueptr.reset(new CppRenyiSignWalk(std_move_wfa(ww.sharedptr.get().create_nonzero_wavefunctionamplitude(ww.sharedptr, deref(rng.uniqueptr.get()))), subsystem.sharedptr))
        return walk

class RenyiSignMeasurementPlan(BasicMeasurementPlan):
    __slots__ = ("walk_plan",)

    def __init__(self, wavefunction, subsystem):
        walk_plan = RenyiSignWalkPlan(wavefunction, subsystem)
        super(RenyiSignMeasurementPlan, self).__init__(walk_plan)

    def to_json(self):
        return OrderedDict([
            ("type", self.__class__.__name__),
            ("subsystem", self.walk_plan.subsystem.to_json()),
        ])

    @staticmethod
    def _from_json(json_repr, wavefunction):
        assert json_repr["type"] == "RenyiSignMeasurementPlan"
        return RenyiSignMeasurementPlan(wavefunction, Subsystem.from_json(json_repr["subsystem"], wavefunction.lattice))

    def to_measurement(self):
        return RenyiSignMeasurement()

cdef class RenyiSignMeasurement(BaseMeasurement):
    def __init__(self):
        self.sharedptr.reset(new CppRenyiSignMeasurement())

    def get_estimate(self, key=None):
        if key is not None:
            raise KeyError
        return Estimate_from_CppComplexBlockedEstimate((<CppRenyiSignMeasurement*>self.sharedptr.get()).get_estimate())

    def get_estimates(self):
        return {None: self.get_estimate()}

class RenyiEntropyMeasurementPlan(CompositeMeasurementPlan):
    def __init__(self, wavefunction, subsystem, steps_per_measurement=None):
        # steps_per_measurement only applies to the
        # SubsystemOccupationProbabilityMeasurementPlan, since the Renyi
        # measurements must be done on every step
        self.sign_plan = RenyiSignMeasurementPlan(wavefunction, subsystem)
        self.modpossible_plan = RenyiModPossibleMeasurementPlan(wavefunction, subsystem)
        extra_args = (steps_per_measurement,) if (steps_per_measurement is not None) else ()
        self.sop_plan = SubsystemOccupationProbabilityMeasurementPlan(wavefunction, subsystem, *extra_args)

    def get_measurement_plans(self):
        return {self.sign_plan, self.modpossible_plan, self.sop_plan}

    def calculate(self, f, key="s2"):
        from math import log
        return {
            "s2_possible": lambda: -log(self.sop_plan.calculate(f, 2)),
            "s2_mod/possible": lambda: -log(f(self.modpossible_plan)),
            "s2_mod": lambda: -log(f(self.modpossible_plan) * self.sop_plan.calculate(f, 2)),
            "s2_sign": lambda: -log(f(self.sign_plan).real),
            "s2": lambda: -log(f(self.modpossible_plan) * f(self.sign_plan).real * self.sop_plan.calculate(f, 2)),
        }[key]()
