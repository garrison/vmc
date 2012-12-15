import abc
import numbers
import collections
from itertools import chain

from pyvmc.utils.immutable import Immutable
from pyvmc.core.lattice import Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.utils import add_hc, ensure_real

class SiteHop(Immutable):
    __slots__ = ("source", "destination", "species")

    def init_validate(self, source, destination, species):
        assert isinstance(source, LatticeSite)
        assert isinstance(destination, LatticeSite)
        assert isinstance(species, numbers.Integral) and species >= 0
        return (source, destination, species)

    def is_valid_for(self, wavefunction):
        assert isinstance(wavefunction, Wavefunction)
        return (self.source in wavefunction.lattice and
                self.destination in wavefunction.lattice and
                self.species < wavefunction.N_species)

    def to_json(self):
        return collections.OrderedDict([
            ("source", self.source.to_json()),
            ("destination", self.destination.to_json()),
            ("species", self.species),
        ])

class Operator(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_basic_operators(self):
        raise NotImplementedError

    @abc.abstractmethod
    def evaluate(self, context):
        """context is a nested dict where a BasicOperator is the key and the value is the expectation value of that operator"""
        raise NotImplementedError

class BasicOperator(Immutable):
    """A BasicOperator represents anything that can be represented a SiteHop's"""

    __slots__ = ("hops", "sum", "boundary_conditions")

    def init_validate(self, hops, sum, boundary_conditions):
        assert isinstance(hops, collections.Sequence)
        hops = tuple(sorted(hops, key=lambda hop: (hop.species, hop.source, hop.destination)))
        assert all([isinstance(hop, SiteHop) for hop in hops])
        assert isinstance(sum, bool)
        if not sum and boundary_conditions is not None:
            from warnings import warn
            warn("boundary_conditions are needlessly given, as sum==False", RuntimeWarning)
            boundary_conditions = None
        if boundary_conditions is not None:
            boundary_conditions = tuple(boundary_conditions)
            assert valid_boundary_conditions(boundary_conditions, len(boundary_conditions))
        return (hops, sum, boundary_conditions)

    def to_json(self):
        return collections.OrderedDict([
            ("type", self.__class__.__name__),
            ("hops", [hop.to_json() for hop in self.hops]),
            ("sum", self.sum),
            ("boundary_conditions", self.boundary_conditions),
        ])

    # these next two methods exist only so we can pass around Operator's
    # without caring whether they are BasicOperator's or CompositeOperator's
    # (see the Operator abstact base class, which declares these methods).

    def get_basic_operators(self):
        return {self}

    def evaluate(self, context):
        def _evaluate():
            return context[self]
        return _evaluate

Operator.register(BasicOperator)

class CompositeOperator(Immutable):
    """
    A CompositeOperator represents anything that is built from BasicOperators
    (e.g. a linear combination of them)
    """

    __slots__ = ("operators",)

    parameters = ()

    def get_basic_operators(self):
        return set(chain.from_iterable([o.get_basic_operators() for o in self.operators]))

Operator.register(CompositeOperator)

# FIXME: move below things to pyvmc.operators or pyvmc.library.operators

class DensityDensityOperator(CompositeOperator):
    __slots__ = ("operators",)

    def init_validate(self, site1, site2, N_species, sum, boundary_conditions):
        # it would be nice if we didn't have to specify N_species here (or at
        # least if we had a way of asserting it's correct in evaluate())...
        assert isinstance(site1, LatticeSite)
        assert isinstance(site2, LatticeSite)
        operators = tuple([BasicOperator([SiteHop(site1, site1, i)] + ([SiteHop(site2, site2, j)] if site1 != site2 or i != j else []), sum, boundary_conditions) for i in xrange(N_species) for j in xrange(N_species)])
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            return sum([o.evaluate(context)() for o in self.operators])
        return _evaluate

class SpinSpinOperator(CompositeOperator):
    __slots__ = ("operators", "samesite")

    def init_validate(self, site1, site2, sum, boundary_conditions):
        if site1 == site2:
            operators = (
                BasicOperator([SiteHop(site1, site1, 0)], sum, boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 1)], sum, boundary_conditions),
            )
            return (operators, True)
        else:
            operators = (
                BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 0)], sum, boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 1)], sum, boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 0), SiteHop(site2, site2, 1)], sum, boundary_conditions),
                BasicOperator([SiteHop(site1, site1, 1), SiteHop(site2, site2, 0)], sum, boundary_conditions),
                BasicOperator([SiteHop(site1, site2, 0), SiteHop(site2, site1, 1)], sum, boundary_conditions),
            )
            return (operators, False)

    def evaluate(self, context):
        operators = self.operators
        def _evaluate():
            if self.samesite:
                return .75 * (context[operators[0]] + context[operators[1]])  # i.e. .75 * rho, scaled by number of sites if sum==True
            else:
                return ensure_real(
                    -.5 * add_hc(context[operators[4]]) +
                    +.25 * context[operators[0]] +
                    +.25 * context[operators[1]] +
                    -.25 * context[operators[2]] +
                    -.25 * context[operators[3]]
                )
        return _evaluate

# FIXME: move everything below to pyvmc.library.dmetal

class SingletRingExchangeOperator(CompositeOperator):
    __slots__ = ("operators",)

    def init_validate(self, sum, boundary_conditions):
        operators = (
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 1)], sum, boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 0)], sum, boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 1)], sum, boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 0)], sum, boundary_conditions),
        )
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            return .5 * sum([o.evaluate(context)() for o in self.operators])
        return _evaluate

class TJKHamiltonian(CompositeOperator):
    __slots__ = ("operators", "divx", "divy")

    parameters = ('t', 'tperp', 'J', 'Jperp', 'K', 'Kperp')

    def init_validate(self, boundary_conditions, lattice):
        assert isinstance(lattice, Lattice)
        assert valid_boundary_conditions(boundary_conditions, len(lattice.dimensions))
        operators = (
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0)], True, boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 1)], True, boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0)], True, boundary_conditions),
            BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 1)], True, boundary_conditions),
            SpinSpinOperator(LatticeSite((0, 0)), LatticeSite((1, 0)), True, boundary_conditions),
            SpinSpinOperator(LatticeSite((0, 0)), LatticeSite((0, 1)), True, boundary_conditions),
            SingletRingExchangeOperator(True, boundary_conditions),
        )

        # don't double-count plaquettes on a 2-leg ladder
        # (FIXME: fix this logic if we migrate to allow cylindrical boundary conditions)
        divx = 2 if (lattice.dimensions[0] == 2 and boundary_conditions is not None) else 1
        divy = 2 if (lattice.dimensions[1] == 2 and boundary_conditions is not None) else 1

        return (operators, divx, divy)

    def evaluate(self, context):
        operators = self.operators
        greenx = context[operators[0]] + context[operators[1]]
        greeny = context[operators[2]] + context[operators[3]]
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

            return (
                -t * add_hc(greenx) / divx +
                -tperp * add_hc(greeny) / divy +
                J * spinspinx / divx +
                Jperp * spinspiny / divy +
                2 * K * add_hc(ringexchange) / divx / divy
            )

        return _evaluate

# FIXME: move below things somewhere else

def _enumerate_with_parity(values):
    """This gives us the correct series of signs if we e.g. multiply out (a-b)(c-d)...(y-z)
    >>> list(_enumerate_with_parity(range(16)))
    [(1, 0), (-1, 1), (-1, 2), (1, 3), (-1, 4), (1, 5), (1, 6), (-1, 7), (-1, 8), (1, 9), (1, 10), (-1, 11), (1, 12), (-1, 13), (-1, 14), (1, 15)]
    """
    for n, value in enumerate(values):
        hamming_weight = bin(n).count('1')
        yield (1 if (hamming_weight % 2 == 0) else -1), value

class SpinSpinTimesSpinSpinOperator(CompositeOperator):
    """Calculates $(S_1 \cdot S_2) (S_3 \cdot S_4)$"""

    __slots__ = ("operators",)

    def init_validate(self, site1, site2, site3, site4, sum, boundary_conditions):
        operators = tuple(chain(
            [
                BasicOperator([
                    SiteHop(site1, site1, i),
                    SiteHop(site2, site2, j),
                    SiteHop(site3, site3, k),
                    SiteHop(site4, site4, l),
                ], sum, boundary_conditions)
                for i in range(2)
                for j in range(2)
                for k in range(2)
                for l in range(2)
            ],
            [
                BasicOperator([
                    SiteHop(site1, site2, 0),
                    SiteHop(site2, site1, 1),
                    SiteHop(site3, site3, i),
                    SiteHop(site4, site4, j),
                ], sum, boundary_conditions)
                for i in range(2)
                for j in range(2)
            ],
            [
                BasicOperator([
                    SiteHop(site3, site4, 0),
                    SiteHop(site4, site3, 1),
                    SiteHop(site1, site1, i),
                    SiteHop(site2, site2, j),
                ], sum, boundary_conditions)
                for i in range(2)
                for j in range(2)
            ],
            [
                BasicOperator([
                    SiteHop(site1, site2, 0),
                    SiteHop(site2, site1, 1),
                    SiteHop(site3, site4, 0),
                    SiteHop(site4, site3, 1),
                ], sum, boundary_conditions),
                BasicOperator([
                    SiteHop(site1, site2, 0),
                    SiteHop(site2, site1, 1),
                    SiteHop(site3, site4, 1),
                    SiteHop(site4, site3, 0),
                ], sum, boundary_conditions),
            ]
        ))
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            assert len(self.operators) == 26
            return sum(chain(
                [.0625 * p * ensure_real(o.evaluate(context)()) for p, o in _enumerate_with_parity(self.operators[0:16])],
                [-.125 * p * add_hc(o.evaluate(context)()) for p, o in _enumerate_with_parity(self.operators[16:20])],
                [-.125 * p * add_hc(o.evaluate(context)()) for p, o in _enumerate_with_parity(self.operators[20:24])],
                [.25 * add_hc(o.evaluate(context)()) for o in self.operators[24:26]],
            ))
        return _evaluate

class SpinModelRingExchangeOperator(CompositeOperator):
    __slots__ = ("operators",)

    def init_validate(self, site1, site2, site3, site4, sum, boundary_conditions):
        """sites 1,2,3,4 should be in cyclical order.

        This is a hermitian operator.
        """
        operators = (
            SpinSpinOperator(site1, site2, sum, boundary_conditions),
            SpinSpinOperator(site1, site3, sum, boundary_conditions),
            SpinSpinOperator(site1, site4, sum, boundary_conditions),
            SpinSpinOperator(site2, site3, sum, boundary_conditions),
            SpinSpinOperator(site2, site4, sum, boundary_conditions),
            SpinSpinOperator(site3, site4, sum, boundary_conditions),
            SpinSpinTimesSpinSpinOperator(site1, site2, site3, site4, sum, boundary_conditions),
            SpinSpinTimesSpinSpinOperator(site1, site4, site2, site3, sum, boundary_conditions),
            SpinSpinTimesSpinSpinOperator(site1, site3, site2, site4, sum, boundary_conditions),
        )
        return (operators,)

    def evaluate(self, context):
        def _evaluate():
            return (
                .25 +
                # two site terms
                self.operators[0].evaluate(context)() +
                self.operators[1].evaluate(context)() +
                self.operators[2].evaluate(context)() +
                self.operators[3].evaluate(context)() +
                self.operators[4].evaluate(context)() +
                self.operators[5].evaluate(context)() +
                # four site terms
                4 * self.operators[6].evaluate(context)() +
                4 * self.operators[7].evaluate(context)() +
                -4 * self.operators[8].evaluate(context)()
            )
        return _evaluate
