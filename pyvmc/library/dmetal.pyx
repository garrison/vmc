import collections
import numbers
from fractions import Fraction

from pyvmc.core.wavefunction import Wavefunction
from pyvmc.core.orbitals import Orbitals
from pyvmc.core.orbitals cimport orbitals_to_orbitaldefinitions

class DMetalWavefunction(Wavefunction):
    """DMetal"""

    __slots__ = ("lattice", "d1", "d2", "f_up", "f_dn", "d1_exponent", "d2_exponent", "f_up_exponent", "f_dn_exponent")

    def init_validate(self, lattice, d1, d2, f_up, f_dn, d1_exponent=1.0, d2_exponent=1.0, f_up_exponent=1.0, f_dn_exponent=1.0):
        (lattice,) = super(DMetalWavefunction, self).init_validate(lattice)
        assert isinstance(d1_exponent, numbers.Real)
        assert isinstance(d2_exponent, numbers.Real)
        assert isinstance(f_up_exponent, numbers.Real)
        assert isinstance(f_dn_exponent, numbers.Real)
        return (
            lattice,
            Orbitals.from_description(d1, lattice),
            Orbitals.from_description(d2, lattice),
            Orbitals.from_description(f_up, lattice),
            Orbitals.from_description(f_dn, lattice),
            float(d1_exponent),
            float(d2_exponent),
            float(f_up_exponent),
            float(f_dn_exponent)
        )

    @property
    def rho(self):
        return Fraction(len(self.d1.momentum_sites), len(self.lattice))

    @property
    def N_species(self):
        return 2

    def to_json(self):
        return collections.OrderedDict([
            ('type', self.__class__.__name__),
            ('orbitals-d1', self.d1.to_json()),
            ('orbitals-d2', self.d2.to_json()),
            ('orbitals-f_up', self.f_up.to_json()),
            ('orbitals-f_dn', self.f_dn.to_json()),
            ('exponent-d1', self.d1_exponent),
            ('exponent-d2', self.d2_exponent),
            ('exponent-f_up', self.f_up_exponent),
            ('exponent-f_dn', self.f_dn_exponent),
        ])

    def to_wavefunction(self):
        cdef WavefunctionWrapper rv = WavefunctionWrapper()
        rv.sharedptr.reset(new CppDMetalWavefunction(
            orbitals_to_orbitaldefinitions(self.d1, self.lattice),
            orbitals_to_orbitaldefinitions(self.d2, self.lattice),
            orbitals_to_orbitaldefinitions(self.f_up, self.lattice),
            orbitals_to_orbitaldefinitions(self.f_dn, self.lattice),
            self.d1_exponent, self.d2_exponent,
            self.f_up_exponent, self.f_dn_exponent
        ))
        return rv

from pyvmc.core.measurement import CompositeMeasurementPlan, BasicOperatorMeasurementPlan
from pyvmc.core.operator import SiteHop, BasicOperator
from pyvmc.core.lattice import LatticeSite
from pyvmc.core.boundary_conditions import periodic # fixme: this is being assumed by these measurements ...

class ElectronRingExchangeMeasurementPlan(CompositeMeasurementPlan):
    __slots__ = ("wavefunction", "plans")
    _immutable_slots = ("wavefunction",)

    def init_validate(self, wavefunction):
        object.__setattr__(self, "plans", (
            BasicOperatorMeasurementPlan(wavefunction, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 1)], (periodic, periodic))),
            BasicOperatorMeasurementPlan(wavefunction, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 0)], (periodic, periodic))),
            BasicOperatorMeasurementPlan(wavefunction, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 1)], (periodic, periodic))),
            BasicOperatorMeasurementPlan(wavefunction, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 0)], (periodic, periodic))),
        ))
        return (wavefunction,)

    def get_measurement_plans(self):
        return set(self.plans)

    def get_result(self, universe):
        return .5 * sum([universe[p].get_result() for p in self.plans]) / len(self.wavefunction.lattice)
