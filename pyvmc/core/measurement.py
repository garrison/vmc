import abc
import numbers
import collections

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.subsystem import Subsystem
from pyvmc.core.lattice import LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
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

class SiteHop(Immutable):
    __slots__ = ("source", "destination", "species")

    def init_validate(self, source, destination, species):
        assert isinstance(source, LatticeSite)
        assert isinstance(destination, LatticeSite)
        assert isinstance(species, numbers.Integral) and species >= 0
        return (source, destination, species)

    def is_valid_for(self, wavefunction):
        assert isinstance(wavefunction, Wavefunction)
        # fixme: implement some actual checks here, including looking at
        # wavefunction.n_species() (which does not yet exist)
        return True

class OperatorMeasurementPlan(MeasurementPlan):
    __slots__ = ("walk", "hops", "sum", "boundary_conditions")

    def __init__(self, wavefunction, hops, sum, boundary_conditions):
        walk = StandardWalkPlan(wavefunction)
        assert isinstance(hops, collections.Sequence)
        hops = tuple(hops)
        assert all([isinstance(hop, SiteHop) and hop.is_valid_for(wavefunction) for hop in hops])
        assert isinstance(sum, bool)
        if boundary_conditions is not None:
            assert sum is True
            assert valid_boundary_conditions(boundary_conditions, len(wavefunction.lattice.dimensions))
            boundary_conditions = tuple(boundary_conditions)
        super(OperatorMeasurementPlan, self).__init__(walk, hops, sum, boundary_conditions)

    def to_json(self):
        lattice = self.walk.wavefunction.lattice
        hops_list = [{
            "source": lattice.index(hop.source),
            "destination": lattice.index(hop.destination),
            "species": hop.species,
        } for hop in self.hops]
        return {
            "type": "operator",
            "hops": hops_list,
            "sum": self.sum,
            "boundary-conditions": self.boundary_conditions,
            "steps-per-measurement": 1,
        }

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
