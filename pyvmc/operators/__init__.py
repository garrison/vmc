from pyvmc.core.lattice import Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.operator import SiteHop, BasicOperator, CompositeOperator
from pyvmc.utils import ensure_real

if True:
    def soft_ensure_hermitian(x):
        return x.real
else:
    soft_ensure_hermitian = ensure_real

class DensityDensityOperator(CompositeOperator):
    __slots__ = ("operators",)

    def init_validate(self, site1, site2, N_species, boundary_conditions):
        # it would be nice if we didn't have to specify N_species here (or at
        # least if we had a way of asserting it's correct in evaluate())...
        assert isinstance(site1, LatticeSite)
        assert isinstance(site2, LatticeSite)
        operators = tuple([BasicOperator([SiteHop(site1, site1, i)] + ([SiteHop(site2, site2, j)] if site1 != site2 or i != j else []), boundary_conditions) for i in range(N_species) for j in range(N_species)])
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            return sum([o.evaluate(context)() for o in self.operators])
        return _evaluate

class SpinSpinOperator(CompositeOperator):
    __slots__ = ("operators", "samesite")

    def init_validate(self, site1, site2, boundary_conditions):
        if site1 == site2:
            operators = (
                BasicOperator([SiteHop(site1, site1, 0)], boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 1)], boundary_conditions),
            )
            return (operators, True)
        else:
            operators = (
                BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 0)], boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 1)], boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 1)], boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 0)], boundary_conditions),
                BasicOperator([SiteHop(site1, site2, 0), SiteHop(site2, site1, 1)], boundary_conditions),
                BasicOperator([SiteHop(site1, site2, 1), SiteHop(site2, site1, 0)], boundary_conditions),
            )
            return (operators, False)

    def evaluate(self, context):
        operators = self.operators
        def _evaluate():
            if self.samesite:
                return .75 * (context[operators[0]] + context[operators[1]])  # i.e. .75 * rho, scaled by number of sites if summing
            else:
                return ensure_real(
                    -.5 * soft_ensure_hermitian(context[operators[4]] + context[operators[5]]) +
                    +.25 * context[operators[0]] +
                    +.25 * context[operators[1]] +
                    -.25 * context[operators[2]] +
                    -.25 * context[operators[3]]
                )
        return _evaluate

# FIXME: move everything below to pyvmc.library.dmetal

class SingletRingExchangeOperator(CompositeOperator):
    __slots__ = ("operators",)

    def init_validate(self, boundary_conditions):
        operators = (
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 1)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 1)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 0)], boundary_conditions),
            # hc
            BasicOperator([SiteHop(LatticeSite((1, 0)), LatticeSite((0, 0)), 0), SiteHop(LatticeSite((0, 1)), LatticeSite((1, 1)), 1)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((1, 0)), LatticeSite((0, 0)), 1), SiteHop(LatticeSite((0, 1)), LatticeSite((1, 1)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 1)), LatticeSite((0, 0)), 0), SiteHop(LatticeSite((1, 0)), LatticeSite((1, 1)), 1)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 1)), LatticeSite((0, 0)), 1), SiteHop(LatticeSite((1, 0)), LatticeSite((1, 1)), 0)], boundary_conditions),
        )
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            return .5 * soft_ensure_hermitian(sum([o.evaluate(context)() for o in self.operators]))
        return _evaluate

class TJKHamiltonian(CompositeOperator):
    __slots__ = ("operators", "divx", "divy")

    parameters = ('t', 'tperp', 'J', 'Jperp', 'K')

    def init_validate(self, boundary_conditions, lattice):
        assert isinstance(lattice, Lattice)
        assert valid_boundary_conditions(boundary_conditions, len(lattice.dimensions))
        operators = (
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 1)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 1)], boundary_conditions),
            SpinSpinOperator(LatticeSite((0, 0)), LatticeSite((1, 0)), boundary_conditions),
            SpinSpinOperator(LatticeSite((0, 0)), LatticeSite((0, 1)), boundary_conditions),
            SingletRingExchangeOperator(boundary_conditions),
            # hc
            BasicOperator([SiteHop(LatticeSite((1, 0)), LatticeSite((0, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((1, 0)), LatticeSite((0, 0)), 1)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 1)), LatticeSite((0, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 1)), LatticeSite((0, 0)), 1)], boundary_conditions),
        )

        # don't double-count plaquettes on a 2-leg ladder
        divx = 2 if (lattice.dimensions[0] == 2 and boundary_conditions[0] != 0) else 1
        divy = 2 if (lattice.dimensions[1] == 2 and boundary_conditions[1] != 0) else 1

        return (operators, divx, divy)

    def evaluate(self, context):
        operators = self.operators
        greenx = soft_ensure_hermitian(context[operators[0]] + context[operators[1]] + context[operators[7]] + context[operators[8]])
        greeny = soft_ensure_hermitian(context[operators[2]] + context[operators[3]] + context[operators[9]] + context[operators[10]])
        spinspinx = operators[4].evaluate(context)()
        spinspiny = operators[5].evaluate(context)()
        ringexchange = operators[6].evaluate(context)()

        divx = self.divx
        divy = self.divy

        def _evaluate(t, J, K, tperp=None, Jperp=None):
            # of course, tperp only means what you think it does if the legs
            # are long in the x direction

            if tperp is None:
                tperp = t
            if Jperp is None:
                Jperp = J

            return ensure_real(
                -t * greenx / divx +
                -tperp * greeny / divy +
                J * spinspinx / divx +
                Jperp * spinspiny / divy +
                2 * K * ringexchange / divx / divy
            )

        return _evaluate

# FIXME: move below things somewhere else

class JKHamiltonian(CompositeOperator):
    __slots__ = ("operators", "divx", "divy")

    parameters = ('J', 'Jperp', 'K')

    def init_validate(self, boundary_conditions, lattice):
        assert isinstance(lattice, Lattice)
        assert valid_boundary_conditions(boundary_conditions, len(lattice.dimensions))
        operators = (
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((1, 0)), LatticeSite((0, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 1)), LatticeSite((0, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 0)], boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((1, 0)), LatticeSite((1, 1)), 0), SiteHop(LatticeSite((0, 1)), LatticeSite((0, 0)), 0)], boundary_conditions),
        )

        # don't double-count plaquettes on a 2-leg ladder
        divx = 2 if (lattice.dimensions[0] == 2 and boundary_conditions[0] != 0) else 1
        divy = 2 if (lattice.dimensions[1] == 2 and boundary_conditions[1] != 0) else 1

        return (operators, divx, divy)

    def evaluate(self, context):
        operators = self.operators
        greenx = soft_ensure_hermitian(context[operators[0]] + context[operators[1]])
        greeny = soft_ensure_hermitian(context[operators[2]] + context[operators[3]])
        ringexchange = soft_ensure_hermitian(operators[4].evaluate(context)() + operators[5].evaluate(context)())

        divx = self.divx
        divy = self.divy

        def _evaluate(J, K, Jperp=None):
            # of course, Jperp only means what you think it does if the legs
            # are long in the x direction

            if Jperp is None:
                Jperp = J

            return ensure_real(
                -J * greenx / divx +
                -Jperp * greeny / divy +
                K * ringexchange / divx / divy
             )

        return _evaluate

# FIXME: move below things somewhere else

class SpinModelRingExchangeOperator(CompositeOperator):
    __slots__ = ("operators",)

    def init_validate(self, site1, site2, site3, site4, boundary_conditions):
        """sites 1,2,3,4 should be in cyclical order.

        This is a hermitian operator.
        """
        operators = (
            # all four spins the same (these operators are hermitian)
            BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 0), SiteHop(site3, site3, 0), SiteHop(site4, site4, 0)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 1), SiteHop(site3, site3, 1), SiteHop(site4, site4, 1)], boundary_conditions),
            # two spins up, two spins down (with all spins being changed)
            BasicOperator([SiteHop(site1, site2, 0), SiteHop(site2, site1, 1), SiteHop(site3, site4, 0), SiteHop(site4, site3, 1)], boundary_conditions),
            BasicOperator([SiteHop(site1, site2, 1), SiteHop(site2, site1, 0), SiteHop(site3, site4, 1), SiteHop(site4, site3, 0)], boundary_conditions),
            # hermitian conjugate
            BasicOperator([SiteHop(site2, site1, 0), SiteHop(site1, site2, 1), SiteHop(site4, site3, 0), SiteHop(site3, site4, 1)], boundary_conditions),
            BasicOperator([SiteHop(site2, site1, 1), SiteHop(site1, site2, 0), SiteHop(site4, site3, 1), SiteHop(site3, site4, 0)], boundary_conditions),
            # three spins one way, one the other
            BasicOperator([SiteHop(site4, site4, 0), SiteHop(site1, site1, 0), SiteHop(site2, site3, 0), SiteHop(site3, site2, 1)], boundary_conditions),
            BasicOperator([SiteHop(site4, site4, 1), SiteHop(site1, site1, 1), SiteHop(site2, site3, 1), SiteHop(site3, site2, 0)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 0), SiteHop(site3, site4, 0), SiteHop(site4, site3, 1)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 1), SiteHop(site3, site4, 1), SiteHop(site4, site3, 0)], boundary_conditions),
            BasicOperator([SiteHop(site2, site2, 0), SiteHop(site3, site3, 0), SiteHop(site4, site1, 0), SiteHop(site1, site4, 1)], boundary_conditions),
            BasicOperator([SiteHop(site2, site2, 1), SiteHop(site3, site3, 1), SiteHop(site4, site1, 1), SiteHop(site1, site4, 0)], boundary_conditions),
            BasicOperator([SiteHop(site3, site3, 0), SiteHop(site4, site4, 0), SiteHop(site1, site2, 0), SiteHop(site2, site1, 1)], boundary_conditions),
            BasicOperator([SiteHop(site3, site3, 1), SiteHop(site4, site4, 1), SiteHop(site1, site2, 1), SiteHop(site2, site1, 0)], boundary_conditions),
            # hermitian conjugate
            BasicOperator([SiteHop(site4, site4, 0), SiteHop(site1, site1, 0), SiteHop(site3, site2, 0), SiteHop(site2, site3, 1)], boundary_conditions),
            BasicOperator([SiteHop(site4, site4, 1), SiteHop(site1, site1, 1), SiteHop(site3, site2, 1), SiteHop(site2, site3, 0)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 0), SiteHop(site4, site3, 0), SiteHop(site3, site4, 1)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 1), SiteHop(site4, site3, 1), SiteHop(site3, site4, 0)], boundary_conditions),
            BasicOperator([SiteHop(site2, site2, 0), SiteHop(site3, site3, 0), SiteHop(site1, site4, 0), SiteHop(site4, site1, 1)], boundary_conditions),
            BasicOperator([SiteHop(site2, site2, 1), SiteHop(site3, site3, 1), SiteHop(site1, site4, 1), SiteHop(site4, site1, 0)], boundary_conditions),
            BasicOperator([SiteHop(site3, site3, 0), SiteHop(site4, site4, 0), SiteHop(site2, site1, 0), SiteHop(site1, site2, 1)], boundary_conditions),
            BasicOperator([SiteHop(site3, site3, 1), SiteHop(site4, site4, 1), SiteHop(site2, site1, 1), SiteHop(site1, site2, 0)], boundary_conditions),
            # two spins up, two spins down (with only two spins being changed)
            BasicOperator([SiteHop(site2, site2, 0), SiteHop(site4, site4, 1), SiteHop(site1, site3, 1), SiteHop(site3, site1, 0)], boundary_conditions),
            BasicOperator([SiteHop(site2, site2, 1), SiteHop(site4, site4, 0), SiteHop(site1, site3, 0), SiteHop(site3, site1, 1)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 0), SiteHop(site3, site3, 1), SiteHop(site4, site2, 1), SiteHop(site2, site4, 0)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 1), SiteHop(site3, site3, 0), SiteHop(site4, site2, 0), SiteHop(site2, site4, 1)], boundary_conditions),
            # hermitian conjugate
            BasicOperator([SiteHop(site2, site2, 0), SiteHop(site4, site4, 1), SiteHop(site3, site1, 1), SiteHop(site1, site3, 0)], boundary_conditions),
            BasicOperator([SiteHop(site2, site2, 1), SiteHop(site4, site4, 0), SiteHop(site3, site1, 0), SiteHop(site1, site3, 1)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 0), SiteHop(site3, site3, 1), SiteHop(site2, site4, 1), SiteHop(site4, site2, 0)], boundary_conditions),
            BasicOperator([SiteHop(site1, site1, 1), SiteHop(site3, site3, 0), SiteHop(site2, site4, 0), SiteHop(site4, site2, 1)], boundary_conditions),
        )
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            assert len(self.operators) == 30
            return soft_ensure_hermitian(2 * sum(o.evaluate(context)() for o in self.operators[0:2])
                                         + sum(o.evaluate(context)() for o in self.operators[2:6])
                                         - sum(o.evaluate(context)() for o in self.operators[6:30]))
        return _evaluate
