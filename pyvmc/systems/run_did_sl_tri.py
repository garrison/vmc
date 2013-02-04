#!/usr/bin/env python

from __future__ import division

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic, boundary_condition_to_string
from math import sqrt, pi

import numpy

import logging
logger = logging.getLogger(__name__)

import os


def alpha_scan(tolerance=None):
    from pyvmc.library.mft_utils import did_nn_bcs_theory, plot_pairing_matrix

    from pyvmc.library.bcs import ProjectedBCSWavefunction

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.core.measurement import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.tmp.scan import do_calculate_plans

    datadir = 'data'
    prefix = 'did_nn-test'

    boundary_conditions = (periodic, antiperiodic)  # for partons in mft (physical bound. conds. will always be fully periodic)
    lattice = HexagonalLattice([18, 18])

    alphas = numpy.arange(pi/100, pi/2, pi/100)

    boundary_conditions_string = [boundary_condition_to_string(bc)[0] for bc in boundary_conditions]
    filename = (prefix + '_' +
               'Lx' + str(lattice.dimensions[0]) + '_' +
               'Ly' + str(lattice.dimensions[1]) + '_' +
               'BCX' + boundary_conditions_string[0] + '_' +
               'BCY' + boundary_conditions_string[1] + '.dat')
    datafile = open(os.path.join(datadir, filename), 'a')

    for alpha in alphas:
        logger.info(' Starting new state, alpha = %f', alpha)

        wf_params = {
            't1': numpy.cos(alpha),
            'delta1': numpy.sin(alpha),
            'delta0': 0,
            'mu0_start': 0,
            'mu0': None,
            'norm': len(lattice)
            }

        bcs_theory = did_nn_bcs_theory(lattice, boundary_conditions, **wf_params)

        phi = bcs_theory['pairing_matrix']

        mu0 = bcs_theory['chemical_potential']
        fdagf = numpy.mean(bcs_theory['bcs_stats']['fdagf'])
        Ttot = numpy.mean(bcs_theory['bcs_stats']['Ttot'])

        wf = ProjectedBCSWavefunction(**{
            'lattice': lattice,
            'phi': phi,
            'N_up': len(lattice) // 2,
        })

        hamiltonian = HeisenbergPlusRingExchangeHamiltonian((periodic, periodic), lattice)
        plans = [BasicOperatorMeasurementPlan(wf, o) for o in hamiltonian.get_basic_operators()]
        results = do_calculate_plans(plans)
        # result[-1] gets the last element of the binned array (FIXME: how to do
        # this? That is, should do_calculate_plans return a data stream or a
        # result?)
        context = {p.operator: result[-1] for p, result in results.iteritems()}
        evaluator = hamiltonian.evaluate(context)

        HeisNN = evaluator(J1=1, J2=0, J3=0, K=0) / len(lattice) / 3
        HeisNNN = evaluator(J1=0, J2=1, J3=0, K=0) / len(lattice) / 3
        HeisNNNN = evaluator(J1=0, J2=0, J3=1, K=0) / len(lattice) / 3
        ring4site = evaluator(J1=0, J2=0, J3=0, K=1) / len(lattice) / 3

        logger.info(' HeisNN = %.6f', HeisNN)
        logger.info(' HeisNNN = %.6f', HeisNNN)
        logger.info(' HeisNNNN = %.6f', HeisNNNN)
        logger.info(' ring4site = %.6f', ring4site)

        data2write = numpy.array([alpha, HeisNN, HeisNNN, HeisNNNN, ring4site, mu0, fdagf, Ttot])
        data2write.tofile(datafile, '  ')
        datafile.write('\n')
        datafile.flush()

    datafile.close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    alpha_scan()
