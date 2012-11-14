from pyvmc.core.measurement import BaseMeasurementPlan, OperatorMeasurementPlan, SiteHop
from pyvmc.core.lattice import LatticeSite
from pyvmc.core.boundary_conditions import periodic # fixme: this is being assumed by these measurements ...
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.utils import add_hc

class GreenMeasurementPlan(BaseMeasurementPlan):
    __slots__ = ("wavefunction", "site1", "site2", "species")

    def init_validate(self, wavefunction, site1, site2, species):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(site1, LatticeSite)
        assert isinstance(site2, LatticeSite)
        assert species < wavefunction.N_species
        return (wavefunction, site1, site2, species)

    def get_measurement_plans(self):
        return {OperatorMeasurementPlan(self.wavefunction, [SiteHop(self.site1, self.site2, self.species)], True, (periodic, periodic))}

    def get_result(self, universe):
        return universe[self.get_measurement_plans().pop()].get_result() / len(self.wavefunction.lattice)

class DensityDensityMeasurementPlan(BaseMeasurementPlan):
    __slots__ = ("wavefunction", "site1", "site2", "plans")
    _immutable_slots = ("wavefunction", "site1", "site2")

    def init_validate(self, wavefunction, site1, site2):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(site1, LatticeSite)
        assert isinstance(site2, LatticeSite)
        object.__setattr_(self, "plans", [OperatorMeasurementPlan(wavefunction, [SiteHop(site1, site1, i)] + ([SiteHop(site2, site2, j)] if site1 != site2 or i != j else []), True, (periodic, periodic)) for i in xrange(wavefunction.N_species) for j in xrange(wavefunction.N_species)])
        return (wavefunction, site1, site2)

    def get_measurement_plans(self):
        return set(self.plans)

    def get_result(self, universe):
        return sum([universe[m].get_result() for m in self.plans]) / len(self.wavefunction.lattice)

class SpinSpinMeasurementPlan(BaseMeasurementPlan):
    __slots__ = ("wavefunction", "site1", "site2", "plans")
    _immutable_slots = ("wavefunction", "site1", "site2")

    def init_validate(self, wavefunction, site1, site2):
        assert isinstance(wavefunction, Wavefunction)
        assert isinstance(site1, LatticeSite)
        assert isinstance(site2, LatticeSite)
        assert wavefunction.N_species == 2
        # fixme: also assert that it's SU(2) ...
        if site1 == site2:
            # same-site anti-commutation relations give a different result
            object.__setattr__(self, "plans", ())
        else:
            object.__setattr__(self, "plans", (
                OperatorMeasurementPlan(wavefunction, [SiteHop(site1, site1, 0), SiteHop(site2, site2, 0)], True, (periodic, periodic)),
                OperatorMeasurementPlan(wavefunction, [SiteHop(site1, site1, 1), SiteHop(site2, site2, 1)], True, (periodic, periodic)),
                OperatorMeasurementPlan(wavefunction, [SiteHop(site1, site1, 0), SiteHop(site2, site2, 1)], True, (periodic, periodic)),
                OperatorMeasurementPlan(wavefunction, [SiteHop(site1, site1, 1), SiteHop(site2, site2, 0)], True, (periodic, periodic)),
                OperatorMeasurementPlan(wavefunction, [SiteHop(site1, site2, 0), SiteHop(site2, site1, 1)], True, (periodic, periodic)),
            ))
        return (wavefunction, site1, site2)

    def get_measurement_plans(self):
        return set(self.plans)

    def get_result(self, universe):
        if self.site1 == self.site2:
            # same-site anti-commutation relations give a different result
            return 0.75 * self.wavefunction.rho
        else:
            return (-.5 * add_hc(universe[self.plans[4]].get_result())
                    + .25 * sum([universe[self.plans[i]].get_result() for i in xrange(0, 2)])
                    - .25 * sum([universe[self.plans[i]].get_result() for i in xrange(2, 4)])) / len(self.wavefunction.lattice)
