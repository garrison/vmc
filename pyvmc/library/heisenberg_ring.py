from itertools import chain

from pyvmc.core.lattice import LatticeRealization, LatticeSite
from pyvmc.core.operator import CompositeOperator, SpinSpinOperator, SpinModelRingExchangeOperator

class HeisenbergPlusRingExchangeHamiltonian(CompositeOperator):
    __slots__ = ("operators", "operators_J1", "operators_J2", "operators_J3", "operators_K")

    def init_validate(self, boundary_conditions, lattice):
        assert isinstance(lattice, LatticeRealization)
        assert len(lattice.dimensions) == 2
        assert lattice.basis_indices == 1

        origin = LatticeSite([2, 2])  # FIXME ! ?
        operators_J1 = tuple([SpinSpinOperator(origin, site, boundary_conditions)
                              for site in lattice.nearest_neighbors(origin, double_count=False)])
        operators_J2 = tuple([SpinSpinOperator(origin, site, boundary_conditions)
                              for site in lattice.second_nearest_neighbors(origin, double_count=False)])
        operators_J3 = tuple([SpinSpinOperator(origin, site, boundary_conditions)
                              for site in lattice.third_nearest_neighbors(origin, double_count=False)])
        operators_K = tuple([SpinModelRingExchangeOperator(site1, site2, site3, site4, boundary_conditions)
                             for site1, site2, site3, site4 in lattice.basic_plaquettes()])

        operators = tuple(chain(
            operators_J1,
            operators_J2,
            operators_J3,
            operators_K
        ))
        return (operators, operators_J1, operators_J2, operators_J3, operators_K)

    def evaluate(self, context):
        def _evaluate(J1, J2, J3, K):
            return (
                J1 * sum(o.evaluate(context)() for o in self.operators_J1) +
                J2 * sum(o.evaluate(context)() for o in self.operators_J2) +
                J3 * sum(o.evaluate(context)() for o in self.operators_J3) +
                K * sum(o.evaluate(context)() for o in self.operators_K)
            )
        return _evaluate
