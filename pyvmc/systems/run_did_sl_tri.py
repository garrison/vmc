#!/usr/bin/env python

from __future__ import division

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic, boundary_condition_to_string
from math import sqrt, pi

import numpy

import logging
logger = logging.getLogger(__name__)

import os

from collections import OrderedDict


def parameter_scan(states_iterable):
    from pyvmc.library.mft_utils import did_nn_bcs_theory, plot_pairing_matrix

    from pyvmc.library.bcs import ProjectedBCSWavefunction

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.core.measurement import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.tmp.scan import do_calculate_plans

    datadir = 'data'
    prefix = 'did_nn-test2'

    lattice = HexagonalLattice([12, 12])
    boundary_conditions = (periodic, antiperiodic)  # for partons in mft (physical bound. conds. will always be fully periodic)

    wf_params = {
        'delta0': 0,
        'mu0_start': 0.1,  # FIXME: be careful with our root-finder at the moment; it's sensitive to mu0_start..
    }

    boundary_conditions_string = [boundary_condition_to_string(bc) for bc in boundary_conditions]
    filename = (prefix + '_' +
               'Lx' + str(lattice.dimensions[0]) + '_' +
               'Ly' + str(lattice.dimensions[1]) + '_' +
               'BCX' + boundary_conditions_string[0][0] + '_' +
               'BCY' + boundary_conditions_string[1][0] + '.dat')
    fullpath = os.path.join(datadir, filename)

    info_keys = (
        'HeisNN',
        'HeisNNN',
        'HeisNNNN',
        'ring4site',
        't1',
        'delta1',
        'delta0',
        'mu0',
        'fdagf',
        'Ttot',
    )

    if not os.path.exists(fullpath):
        datafile = open(fullpath, 'w')
        datafile.write('# lattice dimensions = ' + str(lattice.dimensions) + '\n')
        datafile.write('# boundary conditions = ' + str(boundary_conditions_string) + '\n')
        datafile.write('# ')
        for key in info_keys: datafile.write(key + '  ')
        datafile.write('\n')
    else:
        datafile = open(fullpath, 'a')


    for wf_params_override in states_iterable:
        for k, v in wf_params_override.items():
            wf_params[k] = v
        wf_params['norm'] = len(lattice)  # FIXME..

        logger.info('starting a new state! wf_params = %s', wf_params)
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

        Hami_terms = OrderedDict([
            ('HeisNN', evaluator(J1=1, J2=0, J3=0, K=0) / len(lattice) / 3),
            ('HeisNNN', evaluator(J1=0, J2=1, J3=0, K=0) / len(lattice) / 3),
            ('HeisNNNN', evaluator(J1=0, J2=0, J3=1, K=0) / len(lattice) / 3),
            ('ring4site', evaluator(J1=0, J2=0, J3=0, K=1) / len(lattice) / 3),
        ])

        for term in Hami_terms.items():
            logger.info(term[0] + ' = %.6f', term[1])

        output = {
            'HeisNN': Hami_terms['HeisNN'],
            'HeisNNN': Hami_terms['HeisNNN'],
            'HeisNNNN': Hami_terms['HeisNNNN'],
            'ring4site': Hami_terms['ring4site'],
            't1': wf_params['t1'],
            'delta1': wf_params['delta1'],
            'delta0': wf_params['delta0'],
            'mu0': bcs_theory['chemical_potential'],
            'fdagf': bcs_theory['bcs_stats']['fdagf'],
            'Ttot': bcs_theory['bcs_stats']['Ttot'],
        }

        data2write = numpy.array( [ output[k] for k in info_keys ] )
        data2write.tofile(datafile, '  ')
        datafile.write('\n')
        datafile.flush()

#        wf_params['mu0_start'] = bcs_theory['chemical_potential']

    datafile.close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    def nn_alpha_parameters():
        for alpha in numpy.arange(pi/100, pi/2, pi/100):
            logger.info(' Starting new state, alpha = %f', alpha)
            yield {
                't1': numpy.cos(alpha),
                'delta1': numpy.sin(alpha),
            }

    #parameter_scan(nn_alpha_parameters())

    from pyvmc.utils.parameter_iteration import iterate_parameters
    parameter_scan([d for d in iterate_parameters(['t1', 'delta1']) if d['delta1'] != 0])
