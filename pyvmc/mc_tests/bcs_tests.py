#!/usr/bin/env python

import logging
from math import sqrt

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic

logger = logging.getLogger(__name__)

def test_ansatz(lattice, t1, delta1, mu0, expected_results, tolerance):
    parton_boundary_conditions = (periodic, antiperiodic)

    from pyvmc.library.mft_utils import did_hf_bcs_theory
    phi = did_hf_bcs_theory(lattice, parton_boundary_conditions, t1=t1, delta1=delta1, mu0=mu0, delta0=0)['pairing_matrix']

    from pyvmc.library.bcs import ProjectedBCSWavefunction
    wf = ProjectedBCSWavefunction(**{
        'lattice': lattice,
        'phi': phi,
        'N_up': len(lattice) // 2,
    })

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.measurements import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.core.universe import SimulationUniverse

    hamiltonian = HeisenbergPlusRingExchangeHamiltonian(parton_boundary_conditions, lattice)
    plans = [BasicOperatorMeasurementPlan(wf, o, steps_per_measurement=100) for o in hamiltonian.get_basic_operators()]
    universe = SimulationUniverse(plans, equilibrium_sweeps=100000)
    universe.iterate(1000000)
    context = {mp.operator: m.get_estimate().result for mp, m in universe.get_overall_measurement_dict().items()}
    evaluator = hamiltonian.evaluate(context)
    final_results = (
        evaluator(J1=1, J2=0, J3=0, K=0) / len(lattice) / 3,
        evaluator(J1=0, J2=1, J3=0, K=0) / len(lattice) / 3,
        evaluator(J1=0, J2=0, J3=1, K=0) / len(lattice) / 3,
        evaluator(J1=0, J2=0, J3=0, K=1) / len(lattice) / 3,
    )
    logger.info("t1=%f, delta1=%f, mu0=%f", t1, delta1, mu0)
    logger.info("  nn: %f", final_results[0])
    logger.info(" 2nn: %f", final_results[1])
    logger.info(" 3nn: %f", final_results[2])
    logger.info("ring: %f", final_results[3])

    assert all(abs(r1 - r2) < tolerance for r1, r2 in zip(final_results, expected_results))

def test_projected_bcs_states(tolerance=.015):
    # Z2A
    test_ansatz(HexagonalLattice([6, 6]), t1=0, delta1=1, mu0=0,
                expected_results=(-.1589, .0266, .0279, -.09), tolerance=tolerance)
    #test_ansatz(HexagonalLattice([18, 18]), t1=0, delta1=1, mu0=0,
    #            expected_results=(-.157, .028, .025, -.10), tolerance=tolerance)

    # ChSU2B
    test_ansatz(HexagonalLattice([6, 6]), t1=1, delta1=sqrt(2), mu0=0,
                expected_results=(-.1604, .029, .017, .102), tolerance=tolerance)
    #test_ansatz(HexagonalLattice([18, 18]), t1=1, delta1=sqrt(2), mu0=0,
    #            expected_results=(-.160, .029, .017, .102), tolerance=tolerance)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_projected_bcs_states()
