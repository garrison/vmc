from __future__ import division

from pyvmc.core import HexagonalLattice, Bands, periodic, antiperiodic
from pyvmc.utils import average_and_stddevmean
from math import sqrt, pi

import numpy

import logging
logger = logging.getLogger(__name__)

import os

from collections import OrderedDict


def parameter_scan(theory_func, states_iterable, datadir, prefix, lattice, boundary_conditions, wf_params, use_prev_mu=False):
    from pyvmc.library.bcs import ProjectedBCSWavefunction

    from pyvmc.library.heisenberg_ring import HeisenbergPlusRingExchangeHamiltonian
    from pyvmc.measurements import BasicOperatorMeasurementPlan
    from pyvmc.core import LatticeSite
    from pyvmc.core.universe import do_calculate_plans

    boundary_conditions_string = [str(bc) for bc in boundary_conditions]
    fileroot = (prefix + '_' +
               'Lx' + str(lattice.dimensions[0]) + '_' +
               'Ly' + str(lattice.dimensions[1]) + '_' +
               'BCX' + boundary_conditions_string[0][0] + '_' +
               'BCY' + boundary_conditions_string[1][0])

    logpath = os.path.join(datadir, fileroot + '.log')
    hdlr = logging.FileHandler(logpath)
    hdlr.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
    logging.getLogger().addHandler(hdlr)

    datapath = os.path.join(datadir, fileroot + '.dat')

    info_keys = (
        't1',
        'delta1',
        't2',
        'delta2',
        't3',
        'delta3',
        'delta0',
        'mu0',
        'fdagf',
        'Ttot',
        'HeisNN', 'HeisNN_error',
        'HeisNNN', 'HeisNNN_error',
        'HeisNNNN', 'HeisNNNN_error',
        'ring4site', 'ring4site_error',
    )

    if not os.path.exists(datapath):
        datafile = open(datapath, 'w')
        datafile.write('# lattice dimensions = ' + str(lattice.dimensions) + '\n')
        datafile.write('# boundary conditions = ' + str(boundary_conditions_string) + '\n')
        datafile.write('# ')
        for key in info_keys: datafile.write(key + '  ')
        datafile.write('\n')
        datafile.flush()
    else:
        # FIXME: this should check that boundary conditions and lattice dimensions in datapath actually match current job
        datafile = open(datapath, 'a')


    for wf_params_override in states_iterable:
        try:
            for k, v in wf_params_override.items():
                wf_params[k] = v

            logger.info('Starting a new state! wf_params = %s', wf_params)
            bcs_theory = theory_func(lattice, boundary_conditions, **wf_params)

            phi = bcs_theory['pairing_matrix']
            wf = ProjectedBCSWavefunction(**{
                'lattice': lattice,
                'phi': phi,
                'N_up': len(lattice) // 2,
            })

            hamiltonian = HeisenbergPlusRingExchangeHamiltonian((periodic, periodic), lattice)
            plans = [BasicOperatorMeasurementPlan(wf, o, steps_per_measurement=200) for o in hamiltonian.get_basic_operators()]
            results = do_calculate_plans(plans, equilibrium_sweeps=100000, bins=30, measurement_sweeps_per_bin=100000)
            result_lengths = [len(result) for result in results.itervalues()]
            assert len(set(result_lengths)) == 1
            evaluators = []
            for i in xrange(result_lengths[0]):
                context = {p.operator: result[i] for p, result in results.iteritems()}
                evaluators.append(hamiltonian.evaluate(context))

            # get exchange expectation values *per site*
            Hami_terms = OrderedDict([
                ('HeisNN', average_and_stddevmean([evaluator(J1=1, J2=0, J3=0, K=0) / len(lattice) for evaluator in evaluators])),
                ('HeisNNN', average_and_stddevmean([evaluator(J1=0, J2=1, J3=0, K=0) / len(lattice) for evaluator in evaluators])),
                ('HeisNNNN', average_and_stddevmean([evaluator(J1=0, J2=0, J3=1, K=0) / len(lattice) for evaluator in evaluators])),
                ('ring4site', average_and_stddevmean([evaluator(J1=0, J2=0, J3=0, K=1) / len(lattice) for evaluator in evaluators])),
            ])

            for k, v in Hami_terms.items():
                logger.info('%s = %.6f (%.6f)', k, v[0], v[1])

            # FIXME: also record acceptance rate in data file
            output = {
                't1': wf_params['t1'],
                'delta1': wf_params['delta1'],
                't2': wf_params['t2'],
                'delta2': wf_params['delta2'],
                't3': wf_params['t3'],
                'delta3': wf_params['delta3'],
                'delta0': wf_params['delta0'],
                'mu0': bcs_theory['chemical_potential'],
                'fdagf': bcs_theory['bcs_stats']['fdagf'],
                'Ttot': bcs_theory['bcs_stats']['Ttot'],
                'HeisNN': Hami_terms['HeisNN'][0], 'HeisNN_error': Hami_terms['HeisNN'][1],
                'HeisNNN': Hami_terms['HeisNNN'][0], 'HeisNNN_error': Hami_terms['HeisNNN'][1],
                'HeisNNNN': Hami_terms['HeisNNNN'][0], 'HeisNNNN_error': Hami_terms['HeisNNNN'][1],
                'ring4site': Hami_terms['ring4site'][0], 'ring4site_error': Hami_terms['ring4site'][1],
            }

            data2write = numpy.array( [ output[k] for k in info_keys ] )
            data2write.tofile(datafile, '  ')
            datafile.write('\n')
            datafile.flush()

            if use_prev_mu:
                wf_params['mu0_start'] = bcs_theory['chemical_potential']

        except Exception as e:
            logger.error('Exception raised: %r %r', wf_params, e)

    datafile.close()


def d_nn_alpha_parameters():
    for alpha in numpy.linspace(pi/100, pi/2, 50):
        logger.info('Starting new state, alpha = %f', alpha)
        yield {
            't1': numpy.cos(alpha),
            'delta1': numpy.sin(alpha),
            't2': 0,
            'delta2': 0,
            't3': 0,
            'delta3': 0,
        }

def Samuel_parameters():
    #        for delta1, mu0 in [(5, -0.1), (10, -0.3), (15, -0.5)]:
    #        for delta1, mu0 in [(2.5, 0.1)]:
    for delta1 in [7.5]:
        logger.info('Starting new state, delta1 = %f', delta1)
        yield {
            't1': 1,
            'delta1': delta1,
            't2': 0,
            'delta2': 0,
            't3': 0,
            'delta3': 0,
            'mu0': None,
        }

def Cenke_nnn_alpha_parameters():
    for alpha in numpy.linspace(-pi/2, pi/2, 101):
        logger.info('Starting new state, alpha = %f', alpha)
        yield {
            't1': 0,
            'delta1': numpy.cos(alpha),
            't2': 0,
            'delta2': numpy.sin(alpha),
            't3': 0,
            'delta3': 0,
        }

def sbm_nn_params():
    logger.info('Starting new state!')
    yield {
        't1': 1,
        't2': 0,
        't3': 0,
        'delta1': 0,
        'delta2': 0,
        'delta3': 0,
        'mu0_start': 0,
        'mu0': None,
        'delta0': 0,
    }



if __name__ == "__main__":
    from pyvmc.library.mft_utils import did_hf_bcs_theory, dx2minusy2_hf_bcs_theory, sbm_bcs_theory

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        )
    logging.getLogger().handlers[0].setLevel(logging.ERROR)

#    # template for did
#    parameter_scan(
#       theory_func = did_hf_bcs_theory,
#       states_iterable = d_nn_alpha_parameters(),
##       datadir = '/Users/mishmash/Documents/Research/UCSB/Cenke_SL/Data/did/2013-02-18',
#       datadir = '.',
#       prefix = 'did_nn',
#       lattice = HexagonalLattice([6, 6]),
#       boundary_conditions = (periodic, antiperiodic),  # for partons in mft (physical bound. conds. will always be fully periodic)
#       wf_params = {
#           'delta0': 0,
#           'mu0_start': 0,  # FIXME: be careful with our root-finder at the moment; it's sensitive to mu0_start..
#       },
#       use_prev_mu = True,
#    )

#    # template for nodald
#    parameter_scan(
#        theory_func = dx2minusy2_hf_bcs_theory,
#        states_iterable = d_nn_alpha_parameters(),
##        datadir = '/Users/mishmash/Documents/Research/UCSB/Cenke_SL/Data/nodald/2013-02-18',
#        datadir = '.',
#        prefix = 'nodald_nn',
#        lattice = HexagonalLattice([6, 6]),
#        boundary_conditions = (periodic, antiperiodic),  # for partons in mft (physical bound. conds. will always be fully periodic)
#        wf_params = {
#            'delta0': 0,
#            'mu0_start': 0,  # FIXME: be careful with our root-finder at the moment; it's sensitive to mu0_start..
#        },
#        use_prev_mu = True,
#    )

#    # template for Cenke QBT did state
#    parameter_scan(
#        theory_func = did_hf_bcs_theory,
#        states_iterable = Cenke_nnn_alpha_parameters(),
##        datadir = '/Users/mishmash/Documents/Research/UCSB/Cenke_SL/Data/did/2013-02-18',
#        datadir = '.',
#        prefix = 'did_nnnCenke',
#        lattice = HexagonalLattice([6, 6]),
#        boundary_conditions = (periodic, antiperiodic),  # for partons in mft (physical bound. conds. will always be fully periodic)
#        wf_params = {
#            'mu0': 0,
#            'delta0': 0,
#            'mu0_start': 0,  # FIXME: be careful with our root-finder at the moment; it's sensitive to mu0_start..
#        },
#        use_prev_mu = False,
#    )

    # template for sbm
    parameter_scan(
        theory_func = sbm_bcs_theory,
        states_iterable = sbm_nn_params(),
        datadir = '/Users/mishmash/Documents/Research/UCSB/Cenke_SL/Data/sbm/2013-02-20',
#        datadir = '.',
        prefix = 'sbm_nn',
        lattice = HexagonalLattice([14, 14]),
        boundary_conditions = (periodic, antiperiodic),  # for partons in mft (physical bound. conds. will always be fully periodic)
        wf_params = {},
        use_prev_mu = False,
    )

#    from pyvmc.utils.parameter_iteration import iterate_parameters
#    parameter_scan([d for d in iterate_parameters(['t1', 'delta1']) if d['delta1'] != 0])
