#!/usr/bin/env python

import logging
from math import sqrt

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic
from pyvmc.utils import average

logger = logging.getLogger(__name__)

def test_ansatz(lattice, chi, eta, a0_3, expected_results, tolerance):
    parton_boundary_conditions = (periodic, antiperiodic)

    from pyvmc.library.mft_utils import did_nn_bcs_theory
    phi = did_nn_bcs_theory(lattice, parton_boundary_conditions, t1=chi, delta1=eta, mu0=a0_3, delta0=0)['pairing_matrix']

    from pyvmc.library.bcs import ProjectedBCSWavefunction
    wf = ProjectedBCSWavefunction(**{
        'lattice': lattice,
        'phi': phi,
        'N_up': len(lattice) // 2,
    })

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.core.measurement import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.tmp.scan import do_calculate_plans

    hamiltonian = HeisenbergPlusRingExchangeHamiltonian(parton_boundary_conditions, lattice)
    plans = [BasicOperatorMeasurementPlan(wf, o) for o in hamiltonian.get_basic_operators()]
    results = do_calculate_plans(plans)
    context = {p.operator: average(result) for p, result in results.iteritems()}
    evaluator = hamiltonian.evaluate(context)
    final_results = (
        evaluator(J1=1, J2=0, J3=0, K=0) / len(lattice) / 3,
        evaluator(J1=0, J2=1, J3=0, K=0) / len(lattice) / 3,
        evaluator(J1=0, J2=0, J3=1, K=0) / len(lattice) / 3,
        evaluator(J1=0, J2=0, J3=0, K=1) / len(lattice) / 3,
    )
    logger.info("chi=%f, eta=%f, a0_3=%f", chi, eta, a0_3)
    logger.info("  nn: %f", final_results[0])
    logger.info(" 2nn: %f", final_results[1])
    logger.info(" 3nn: %f", final_results[2])
    logger.info("ring: %f", final_results[3])

    assert all(abs(r1 - r2) < tolerance for r1, r2 in zip(final_results, expected_results))

def test_projected_bcs_states(tolerance=.01):
    # Z2A
    test_ansatz(HexagonalLattice([6, 6]), chi=0, eta=1, a0_3=0,
                expected_results=(-.157, .028, .025, -.10), tolerance=2*tolerance)
    #test_ansatz(HexagonalLattice([18, 18]), chi=0, eta=1, a0_3=0,
    #            expected_results=(-.157, .028, .025, -.10), tolerance=tolerance)

    # ChSU2B
    test_ansatz(HexagonalLattice([6, 6]), chi=1, eta=sqrt(2), a0_3=0,
                expected_results=(-.160, .029, .017, .102), tolerance=2*tolerance)
    #test_ansatz(HexagonalLattice([18, 18]), chi=1, eta=sqrt(2), a0_3=0,
    #            expected_results=(-.160, .029, .017, .102), tolerance=tolerance)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_projected_bcs_states()
