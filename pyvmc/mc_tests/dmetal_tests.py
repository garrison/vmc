#!/usr/bin/env python

import logging

from pyvmc.core import Lattice, Bands, periodic, antiperiodic
from pyvmc.library.dmetal import DMetalWavefunction
from pyvmc.utils import average

logger = logging.getLogger(__name__)

def test_dmetal_energy(tolerance=None):
    wf = DMetalWavefunction(**{
        'lattice': Lattice([12, 2]),
        'd1': Bands([5, 3], (periodic, periodic)),
        'd2': Bands([8, 0], (antiperiodic, periodic)),
        'f_up': Bands([4, 0], (antiperiodic, periodic)),
        'f_dn': Bands([4, 0], (antiperiodic, periodic)),
        'd1_exponent': 0.7,
        'd2_exponent': -0.4,
    })

    from pyvmc.operators import TJKHamiltonian
    from pyvmc.measurements import BasicOperatorMeasurementPlan
    from pyvmc.tmp.scan import do_calculate_plans

    hamiltonian = TJKHamiltonian((periodic, periodic), wf.lattice)
    plans = [BasicOperatorMeasurementPlan(wf, o) for o in hamiltonian.get_basic_operators()]
    results = do_calculate_plans(plans)
    context = {p.operator: average(result) for p, result in results.iteritems()}
    energy = hamiltonian.evaluate(context)(t=1, J=2, K=2) / len(wf.lattice)

    logger.info("Energy: %f", energy)
    assert -0.8 < energy < -0.74

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_dmetal_energy()
