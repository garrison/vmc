#!/usr/bin/env python

import logging

from pyvmc.core import Lattice, Bands, periodic, antiperiodic

logger = logging.getLogger(__name__)

def test_spinmodel(tolerance=None):
    from pyvmc.library.dmetal import DMetalWavefunction

    wf = DMetalWavefunction(**{
        'lattice': Lattice([12, 2]),
        'd1': Bands([5, 3], (periodic, periodic)),
        'd2': Bands([8, 0], (antiperiodic, periodic)),
        'f_up': Bands([4, 0], (antiperiodic, periodic)),
        'f_dn': Bands([4, 0], (antiperiodic, periodic)),
        'd1_exponent': 0.7,
        'd2_exponent': -0.4,
    })

    from pyvmc.operators import SpinModelRingExchangeOperator
    from pyvmc.measurements import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.core.universe import SimulationUniverse

    spin_operator = SpinModelRingExchangeOperator(LatticeSite([0, 0]), LatticeSite([1, 0]),
                                                  LatticeSite([1, 1]), LatticeSite([0, 1]),
                                                  (periodic, periodic))
    plans = [BasicOperatorMeasurementPlan(wf, o) for o in spin_operator.get_basic_operators()]
    universe = SimulationUniverse(plans, equilibrium_sweeps=500000)
    universe.iterate(1000000)
    context = {mp.operator: m.get_estimate().result for mp, m in universe.get_overall_measurement_dict().items()}
    logger.info("Spin model: %f", spin_operator.evaluate(context)())

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_spinmodel()
