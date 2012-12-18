#!/usr/bin/env python

import logging

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic

logger = logging.getLogger(__name__)

def test_spinmodel(tolerance=None):
    boundary_conditions = (periodic, periodic)
    lattice = HexagonalLattice([18, 18])

    from pyvmc.library.dmetal import DMetalWavefunction
    wf = DMetalWavefunction(**{
        'lattice': lattice,
        'd1': Bands([5, 3], (periodic, periodic)),
        'd2': Bands([8, 0], (antiperiodic, periodic)),
        'f_up': Bands([4, 0], (antiperiodic, periodic)),
        'f_dn': Bands([4, 0], (antiperiodic, periodic)),
        'd1_exponent': 0.7,
        'd2_exponent': -0.4,
    })

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.core.measurement import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.tmp.scan import do_calculate_plans

    hamiltonian = HeisenbergPlusRingExchangeHamiltonian(boundary_conditions, lattice)
    plans = [BasicOperatorMeasurementPlan(wf, o) for o in hamiltonian.get_basic_operators()]
    results = do_calculate_plans(plans)
    # result[-1] gets the last element of the binned array (FIXME: how to do
    # this? That is, should do_calculate_plans return a data stream or a
    # result?)
    context = {p.operator: result[-1] for p, result in results.iteritems()}
    logger.info("Energy: %f", hamiltonian.evaluate(context)(J1=0, J2=0, J3=0, K=1))

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_spinmodel()
