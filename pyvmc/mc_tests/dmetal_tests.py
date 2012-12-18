#!/usr/bin/env python

import logging

from pyvmc.core import Lattice, Bands, periodic, antiperiodic
from pyvmc.library.dmetal import DMetalWavefunction

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

    # old method

    import h5py
    from pyvmc.library.dmetal import TJKEnergetics

    e = TJKEnergetics(wf)
    with h5py.File('/tmp/pyvmc_mc_test.hdf5', 'w', libver='latest') as f:
        e.calculate(f)
        e.load_results(f)
        energy = e.get_energy(J=2, K=2)

    logger.info("Energy (old method): %f", energy)
    assert energy < -0.74 and energy > -0.8

    # new method

    from pyvmc.core.operator import TJKHamiltonian
    from pyvmc.core.measurement import BasicOperatorMeasurementPlan
    from pyvmc.tmp.scan import do_calculate_plans

    hamiltonian = TJKHamiltonian((periodic, periodic), wf.lattice)
    plans = [BasicOperatorMeasurementPlan(wf, o) for o in hamiltonian.get_basic_operators()]
    results = do_calculate_plans(plans)
    # result[-1] gets the last element of the binned array (FIXME: how to do
    # this? That is, should do_calculate_plans return a data stream or a
    # result?)
    context = {p.operator: result[-1] for p, result in results.iteritems()}
    energy = hamiltonian.evaluate(context)(t=1, J=2, K=2) / len(wf.lattice)

    logger.info("Energy (new method): %f", energy)
    assert energy < -0.74 and energy > -0.8

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_dmetal_energy()
