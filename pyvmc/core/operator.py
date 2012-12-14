import abc
import numbers
import collections
import itertools

from pyvmc.utils.immutable import Immutable
from pyvmc.core.lattice import Lattice, LatticeSite
from pyvmc.core.boundary_conditions import valid_boundary_conditions
from pyvmc.core.wavefunction import Wavefunction
from pyvmc.utils import add_hc

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
        return set(itertools.chain.from_iterable([o.get_basic_operators() for o in self.operators]))

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
                return (
                    -.5 * add_hc(context[operators[4]]) +
                    +.25 * context[operators[0]] +
                    +.25 * context[operators[1]] +
                    -.25 * context[operators[2]] +
                    -.25 * context[operators[3]]
                ).real
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
