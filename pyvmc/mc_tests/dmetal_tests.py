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

    from pyvmc.operators import TJKHamiltonian
    from pyvmc.measurements import BasicOperatorMeasurementPlan
    from pyvmc.core.universe import SimulationUniverse, save_universe_to_hdf5

    hamiltonian = TJKHamiltonian((periodic, periodic), wf.lattice)
    plans = [BasicOperatorMeasurementPlan(wf, o) for o in hamiltonian.get_basic_operators()]
    universe = SimulationUniverse(plans)
    universe.iterate(1000000)
    context = {mp.operator: m.get_estimate().result for mp, m in universe.get_overall_measurement_dict().items()}
    energy = hamiltonian.evaluate(context)(t=1, J=2, K=2) / len(wf.lattice)

    logger.info("Energy: %f", energy)
    assert -0.8 < energy < -0.74

    # test hdf5
    import h5py
    filename = '/tmp/hdf5test.hdf5'
    with h5py.File(filename, 'w') as f:
        grp = f.create_group('testgroup')
        save_universe_to_hdf5(universe, grp)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_dmetal_energy()
