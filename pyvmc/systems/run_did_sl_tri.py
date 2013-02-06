#!/usr/bin/env python

from __future__ import division

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic, boundary_condition_to_string
from math import sqrt, pi

import numpy

import logging
logger = logging.getLogger(__name__)

import os

from collections import OrderedDict


def alpha_scan(tolerance=None):
    from pyvmc.library.mft_utils import did_nn_bcs_theory, plot_pairing_matrix

    from pyvmc.library.bcs import ProjectedBCSWavefunction

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.core.measurement import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.tmp.scan import do_calculate_plans

    datadir = 'data'
    prefix = 'did_nn-head'

    lattice = HexagonalLattice([12, 12])
    boundary_conditions = (periodic, antiperiodic)  # for partons in mft (physical bound. conds. will always be fully periodic)

    alphas = numpy.arange(pi/100, pi/2, pi/100)

    wf_params = OrderedDict([
        ('t1', None),
        ('delta1', None),
        ('delta0', 0),
        ('mu0_start', 0),
        ('mu0', None),
        ('norm', None)
    ])

    wf_stats = OrderedDict([
        ('fdagf', None),
        ('Ttot', None)
    ])

    Hami_terms = OrderedDict([
        ('HeisNN', None),
        ('HeisNNN', None),
        ('HeisNNNN', None),
        ('ring4site', None)
    ])

    boundary_conditions_string = [boundary_condition_to_string(bc) for bc in boundary_conditions]
    filename = (prefix + '_' +
               'Lx' + str(lattice.dimensions[0]) + '_' +
               'Ly' + str(lattice.dimensions[1]) + '_' +
               'BCX' + boundary_conditions_string[0][0] + '_' +
               'BCY' + boundary_conditions_string[1][0] + '.dat')
    fullpath = os.path.join(datadir, filename)

    info_keys = OrderedDict(Hami_terms.items() + wf_params.items() + wf_stats.items()).keys()

    if not os.path.exists(fullpath):
        datafile = open(fullpath, 'w')
        datafile.write('# lattice dimensions = ' + str(lattice.dimensions) + '\n')
        datafile.write('# boundary conditions = ' + str(boundary_conditions_string) + '\n')
        datafile.write('# ')
        for key in info_keys: datafile.write(key + '  ')
        datafile.write('\n')
    else:
        datafile = open(fullpath, 'a')


    for alpha in alphas:
        logger.info(' Starting new state, alpha = %f', alpha)

        wf_params['t1'] = numpy.cos(alpha)
        wf_params['delta1'] = numpy.sin(alpha)
        wf_params['norm'] = len(lattice)  # FIXME..

        bcs_theory = did_nn_bcs_theory(lattice, boundary_conditions, **wf_params)

        phi = bcs_theory['pairing_matrix']

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

        Hami_terms['HeisNN'] = evaluator(J1=1, J2=0, J3=0, K=0) / len(lattice) / 3
        Hami_terms['HeisNNN'] = evaluator(J1=0, J2=1, J3=0, K=0) / len(lattice) / 3
        Hami_terms['HeisNNNN'] = evaluator(J1=0, J2=0, J3=1, K=0) / len(lattice) / 3
        Hami_terms['ring4site'] = evaluator(J1=0, J2=0, J3=0, K=1) / len(lattice) / 3

        for term in Hami_terms.items():
            logger.info(term[0] + ' = %.6f', term[1])

        wf_params['mu0'] = bcs_theory['chemical_potential']
        wf_stats['fdagf'] = bcs_theory['bcs_stats']['fdagf']
        wf_stats['Ttot'] = bcs_theory['bcs_stats']['Ttot']

        data2write = numpy.array( OrderedDict(Hami_terms.items() + wf_params.items() + wf_stats.items()).values() )
        data2write.tofile(datafile, '  ')
        datafile.write('\n')
        datafile.flush()

        # fix things up for next iteration..
        wf_params['mu0_start'] = wf_params['mu0']
        wf_params['mu0'] = None

    datafile.close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    alpha_scan()
