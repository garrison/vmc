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

from pyvmc.core.measurement import BaseMeasurementPlan, BasicOperatorMeasurementPlan
from pyvmc.core.operator import SiteHop, BasicOperator
from pyvmc.core.lattice import LatticeSite
from pyvmc.core.boundary_conditions import periodic # fixme: this is being assumed by these measurements ...

class ElectronRingExchangeMeasurementPlan(BaseMeasurementPlan):
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

from pyvmc.tmp.scan import calculate_plans, load_results

class TJKEnergetics(object):
    def __init__(self, wavefunction):
        self.wavefunction = wavefunction
        from pyvmc.tmp.old_measurements import GreenMeasurementPlan, SpinSpinMeasurementPlan

        origin = LatticeSite((0, 0))
        self.plans = {
            'greenx_up': GreenMeasurementPlan(wavefunction, origin, LatticeSite((1, 0)), 0),
            'greeny_up': GreenMeasurementPlan(wavefunction, origin, LatticeSite((0, 1)), 0),
            'greenx_down': GreenMeasurementPlan(wavefunction, origin, LatticeSite((1, 0)), 1),
            'greeny_down': GreenMeasurementPlan(wavefunction, origin, LatticeSite((0, 1)), 1),
            'spinspinx': SpinSpinMeasurementPlan(wavefunction, origin, LatticeSite((1, 0))),
            'spinspiny': SpinSpinMeasurementPlan(wavefunction, origin, LatticeSite((0, 1))),
            'ringexchange': ElectronRingExchangeMeasurementPlan(wavefunction),
        }

    def calculate(self, h5group):
        calculate_plans(self.plans, h5group)

    def load_results(self, h5group):
        self.results = load_results(self.plans, h5group)

    def get_energy(self, J, K):
        results = self.results

        # don't double-count plaquettes on a 2-leg lattice
        dimensions = self.wavefunction.lattice.dimensions
        assert len(dimensions) == 2
        assert dimensions[0] > 1 and dimensions[1] > 1
        divx = 2 if dimensions[0] == 2 else 1
        divy = 2 if dimensions[1] == 2 else 1

        from pyvmc.utils import add_hc

        return sum([
            -add_hc(results['greenx_up']) / divx,
            -add_hc(results['greenx_down']) / divx,
            -add_hc(results['greeny_up']) / divy,
            -add_hc(results['greeny_down']) / divy,
            J * results['spinspinx'] / divx,
            J * results['spinspiny'] / divy,
            2 * K * add_hc(results['ringexchange']) / divx / divy,
        ])
