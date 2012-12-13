from twisted.internet import defer

from pyvmc.core import SimpleSubsystem
from pyvmc.core.measurement import SubsystemOccupationProbabilityMeasurementPlan, BasicOperatorMeasurementPlan
from pyvmc.core.operator import BasicOperator, SiteHop
from pyvmc.library.renyi import RenyiSignMeasurementPlan, RenyiModPossibleMeasurementPlan
from pyvmc.tmp.universe import TmpMeasurementPlan, get_default_universe

# there are two types of goals: ones which are measured directly in the VMC,
# and ones which depend on other goals

class Goal(object):
    def __init__(self, independent, universe):
        self.n_independent = independent
        self.universe = universe

class SubsystemOccupationNumberProbability(Goal):
    def __init__(self, system, subsystem, independent=30, universe=None):
        if universe is None:
            universe = get_default_universe()

        measurement = SubsystemOccupationProbabilityMeasurementPlan(system, subsystem)
        plan = TmpMeasurementPlan(measurement, self.parse_result)
        self.measurement_set = universe.get_measurement_set(plan, independent)

        self.subsystem = subsystem

    def advance(self):
        return self.measurement_set.advance()

    @staticmethod
    def parse_result(result_json):
        rv = {}
        for occupation, probability in result_json:
            occupation = tuple(occupation)
            rv[occupation] = probability
        return rv

    def get_swap_probability(self, index, n_swaps=2):
        return sum(v ** n_swaps for v in self.measurement_set[index].get_aggregate_result().values())

class RenyiModPossible(Goal):
    def __init__(self, system, subsystem, independent=30, universe=None):
        if universe is None:
            universe = get_default_universe()

        measurement = RenyiModPossibleMeasurementPlan(system, subsystem)
        plan = TmpMeasurementPlan(measurement, self.parse_result)
        self.measurement_set = universe.get_measurement_set(plan, independent)

        self.subsystem = subsystem

    def advance(self):
        return self.measurement_set.advance()

    @staticmethod
    def parse_result(result_json):
        assert isinstance(result_json, float)
        return result_json

    def get_renyi_mod_possible(self, index):
        from math import log
        return -log(self.measurement_set[index].get_aggregate_result())

class RenyiSign(Goal):
    def __init__(self, system, subsystem, independent=30, universe=None):
        if universe is None:
            universe = get_default_universe()

        measurement = RenyiSignMeasurementPlan(system, subsystem)
        plan = TmpMeasurementPlan(measurement, self.parse_result)
        self.measurement_set = universe.get_measurement_set(plan, independent)

        self.subsystem = subsystem

    def advance(self):
        return self.measurement_set.advance()

    @staticmethod
    def parse_result(result_json):
        assert isinstance(result_json, complex)
        return result_json

    def get_renyi_sign(self, index):
        from math import log
        try:
            return -log(self.measurement_set[index].get_aggregate_result().real)
        except Exception:
            import logging
            # if this turns up negative, then the orbitals probably aren't
            # normalized correct (or some over/under flow is happening)
            logging.error("domain error %s", self.measurement_set[index].get_aggregate_result().real)
            return 0.0

class RenyiMod(Goal):
    def __init__(self, system, subsystem, independent=30, universe=None):
        self.renyi_mod_possible = RenyiModPossible(system, subsystem, independent, universe)
        self.sonp = SubsystemOccupationNumberProbability(system, subsystem, independent, universe)
        self.subsystem = subsystem

    def advance(self):
        d1 = self.renyi_mod_possible.advance()
        d2 = self.sonp.advance()
        d = defer.gatherResults([d1, d2], consumeErrors=True)
        return d

    def get_renyi_mod(self, index):
        return self.get_renyi_mod_possible(index) + self.get_renyi_swap_possible(index)

    def get_renyi_swap_possible(self, index):
        from math import log
        # fixme: could be taking log of zero ...
        return -log(self.sonp.get_swap_probability(index))

    def get_renyi_mod_possible(self, index):
        return self.renyi_mod_possible.get_renyi_mod_possible(index)

class Renyi(Goal):
    """given a system and a subsystem size, we want the renyi entropy of it"""

    def __init__(self, system, subsystem, independent=30, universe=None):
        self.renyi_sign = RenyiSign(system, subsystem, independent, universe)
        self.renyi_mod = RenyiMod(system, subsystem, independent, universe)
        self.subsystem = subsystem

    def advance(self):
        """Returns a deferred"""
        # fixme: handle correctly if an advancement is already queued
        d1 = self.renyi_mod.advance()
        d2 = self.renyi_sign.advance()
        d = defer.gatherResults([d1, d2], consumeErrors=True)
        # we might want to add a callback that ... tells us to clear the cache :)
        return d

    def get_renyi_sign(self, index):
        return self.renyi_sign.get_renyi_sign(index)

    def get_renyi_mod(self, index):
        return self.renyi_mod.get_renyi_mod(index)

    def get_renyi(self, index):
        return self.get_renyi_sign(index) + self.get_renyi_mod(index)

class RenyiLengthScaling(Goal):
    """given a quasi-1d system, we want to figure out c based on the renyi entropies"""

    def __init__(self, system, renyi_class=Renyi, independent=30, universe=None):
        subsystem_lengths = xrange(1, (system.lattice.dimensions[0] + 1) // 2 + 1)
        subsystems = [SimpleSubsystem((x,) + system.lattice.dimensions[1:], system.lattice)
                      for x in subsystem_lengths]
        self.renyi_plans = [renyi_class(system, subsystem, independent, universe)
                            for subsystem in subsystems]
        self.system = system

    def advance(self):
        """Returns a deferred"""
        # fixme: handle correctly if an advancement is already queued
        deferreds = [r.advance() for r in self.renyi_plans]
        assert len(deferreds) > 0  # i'm not sure if it works if there are zero...
        return defer.gatherResults(deferreds, consumeErrors=True)

#    def _get_c(self):
#        class A(object):
#            pass
#        rv = A()
#        rv.value = '373'
#        return rv
#
#    c = property(_get_c, doc="conformal charge c")

    def get_c(self, index):
        system_length = self.system.lattice.dimensions[0]
        def conformal_length(x):
            from math import pi, log, sin
            return log(system_length / pi * sin(pi * x / system_length))
        conformal_lengths = [conformal_length(p.subsystem.dimensions[0]) for p in self.renyi_plans]
        renyi_entropies = [p.get_renyi(index) for p in self.renyi_plans]

        import numpy
        fit = numpy.polyfit(conformal_lengths, renyi_entropies, 1)
        return 4.0 * fit[0]

class DensityDensity(Goal):
    @staticmethod
    def parse_result(result_json):
        rv = {}
        assert isinstance(result_json, list)
        for basis_index, bi_result in enumerate(result_json):
            assert isinstance(bi_result, list)
            for position_index, value in enumerate(bi_result):
                assert isinstance(value, float)
                rv[Point(basis_index, position_index)] = value
        return rv

class Green(Goal):
    @staticmethod
    def parse_result(result_json):
        rv = {}
        assert isinstance(result_json, list)
        for basis_index, bi_result in enumerate(result_json):
            assert isinstance(bi_result, list)
            for position_index, value in enumerate(bi_result):
                rv[Point(basis_index, position_index)] = value
        return rv

class Operator(Goal):
    def __init__(self, system, hops, sum, boundary_conditions, independent=30, universe=None):
        if universe is None:
            universe = get_default_universe()

        measurement = BasicOperatorMeasurementPlan(system, BasicOperator(hops, sum, boundary_conditions))
        plan = TmpMeasurementPlan(measurement, self.parse_result)
        self.measurement_set = universe.get_measurement_set(plan, independent)

    def advance(self):
        return self.measurement_set.advance()

    @staticmethod
    def parse_result(result_json):
        assert isinstance(result_json, complex)
        return result_json

    def get_expectation_value(self, index):
        return self.measurement_set[index].get_aggregate_result()

class SingletPairOperator(Goal):
    def __init__(self, system, r1, r2, r1p, r2p, sum, boundary_conditions, independent=30, universe=None):
        """r1, r2, r1p, r2p should all be unique sites"""
        hops1 = (
            SiteHop(r2p, r2, 1),
            SiteHop(r1p, r1, 0),
        )
        hops2 = (
            SiteHop(r1p, r2, 1),
            SiteHop(r2p, r1, 0),
        )
        hops3 = (
            SiteHop(r2p, r1, 1),
            SiteHop(r1p, r2, 0),
        )
        hops4 = (
            SiteHop(r1p, r1, 1),
            SiteHop(r2p, r2, 0),
        )
        self.operators = (
            Operator(system, hops1, sum, boundary_conditions, independent, universe),
            Operator(system, hops2, sum, boundary_conditions, independent, universe),
            Operator(system, hops3, sum, boundary_conditions, independent, universe),
            Operator(system, hops4, sum, boundary_conditions, independent, universe),
        )

    def advance(self):
        d = defer.gatherResults([o.advance() for o in self.operators], consumeErrors=True)
        return d

    def get_expectation_value(self, index):
        return sum([o.get_expectation_value(index) for o in self.operators])

class DiagonalDWaveCooperPairOperator(Goal):
    def __init__(self, system, x, sum, boundary_conditions, independent=30, universe=None):
        # assert two leg ladder
        assert len(system.lattice.dimensions) == 2 and system.lattice.dimensions[1] == 2
        from pyvmc.core.lattice import LatticeSite
        r1 = LatticeSite([0, 0])
        r2 = LatticeSite([1, 1])
        r1p = LatticeSite([x, 0])
        r2p = LatticeSite([x + 1, 1])
        r1ps = LatticeSite([x, 1])
        r2ps = LatticeSite([x + 1, 0])
        self.operators = (
            SingletPairOperator(system, r1, r2, r1p, r2p, sum, boundary_conditions, independent, universe),
            SingletPairOperator(system, r1, r2, r1ps, r2ps, sum, boundary_conditions, independent, universe),
        )

    def advance(self):
        d = defer.gatherResults([o.advance() for o in self.operators], consumeErrors=True)
        return d

    def get_symmetric_expectation_value(self, index):
        return 2 * (self.operators[0].get_expectation_value(index)
                    + self.operators[1].get_expectation_value(index))

    def get_antisymmetric_expectation_value(self, index):
        return 2 * (self.operators[0].get_expectation_value(index)
                    - self.operators[1].get_expectation_value(index))

    # temporary, for backwards compatibility
    get_expectation_value = get_antisymmetric_expectation_value

class LegDWaveCooperPairOperator(Goal):
    def __init__(self, system, x, sum, boundary_conditions, independent=30, universe=None):
        # assert two leg ladder
        assert len(system.lattice.dimensions) == 2 and system.lattice.dimensions[1] == 2
        from pyvmc.core.lattice import LatticeSite
        r1 = LatticeSite([0, 0])
        r2 = LatticeSite([1, 0])
        r1p = LatticeSite([x, 0])
        r2p = LatticeSite([x + 1, 0])
        r1ps = LatticeSite([x, 1])
        r2ps = LatticeSite([x + 1, 1])
        self.operators = (
            SingletPairOperator(system, r1, r2, r1p, r2p, sum, boundary_conditions, independent, universe),
            SingletPairOperator(system, r1, r2, r1ps, r2ps, sum, boundary_conditions, independent, universe),
        )

    def advance(self):
        d = defer.gatherResults([o.advance() for o in self.operators], consumeErrors=True)
        return d

    def get_symmetric_expectation_value(self, index):
        return 2 * (self.operators[0].get_expectation_value(index)
                    + self.operators[1].get_expectation_value(index))

    def get_antisymmetric_expectation_value(self, index):
        return 2 * (self.operators[0].get_expectation_value(index)
                    - self.operators[1].get_expectation_value(index))

    # temporary, for backwards compatibility
    get_expectation_value = get_antisymmetric_expectation_value

class DensityDensityFourier(Goal):
    pass

class GreenFourier(Goal):
    pass
