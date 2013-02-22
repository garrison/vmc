import six

from itertools import chain
from collections import Sequence
import logging

import numpy

from pyvmc.utils import custom_json as json
from pyvmc.core.simulation import MetropolisSimulation
from pyvmc.core.measurement import MeasurementPlan, BasicMeasurementPlan
from pyvmc.core.rng import RandomNumberGenerator
from pyvmc.utils import average

logger = logging.getLogger(__name__)

def _create_basic_measurement_plan_set(measurement_plans):
    return set(chain.from_iterable(p.get_measurement_plans() for p in measurement_plans))

# fixme: can we guarantee that p.to_json() doesn't have any forward slashes?
# fixme: we should really save each measurement in a subgroup for that walk, and in there store information from the walk's completion... INCLUDING rusage, etc

class SimulationUniverse(object):
    """
    """

    def __init__(self, measurement_plans, equilibrium_sweeps=500000):
        # NOTE: it may one day be useful to keep a list of all
        # measurement_plans originally passed, but at the moment this does not
        # seem to be useful.

        # first create all the basic measurements that need to be performed
        assert isinstance(measurement_plans, Sequence)
        assert all(isinstance(mp, MeasurementPlan) for mp in measurement_plans)
        basic_measurement_plans = list(_create_basic_measurement_plan_set(measurement_plans))
        assert all(isinstance(mp, BasicMeasurementPlan) for mp in basic_measurement_plans)

        # now organize them by each walk which must be performed
        measurement_plans_by_walk = {}
        for mp in basic_measurement_plans:
            measurement_plans_by_walk.setdefault(mp.walk_plan, []).append(mp)

        # prepare and equilibriate simulations
        self.simulations = []
        for walk_plan, current_measurement_plans in six.iteritems(measurement_plans_by_walk):
            # the MetropolisSimulation constructor eats the RNG, so we need to
            # create a new one for each simulation
            rng = RandomNumberGenerator()
            self.simulations.append(MetropolisSimulation(walk_plan, walk_plan.wavefunction.lattice,
                                                         current_measurement_plans, equilibrium_sweeps, rng))

    def iterate(self, sweeps):
        for sim in self.simulations:
            sim.reset_measurement_estimates()
            sim.iterate(sweeps)

        # fixme: should we return anything at all?
        from itertools import chain
        return dict(chain.from_iterable(six.iteritems(sim.measurement_dict) for sim in self.simulations))

def do_calculate_plans(plans, equilibrium_sweeps=500000, bins=100, measurement_sweeps_per_bin=10000):
    # prepare and equilibrate simulations
    calc = SimulationUniverse(plans, equilibrium_sweeps)

    # perform simulations
    results = {}
    for i in six.moves.xrange(bins):  # FIXME: this outer loop should be unnecessary.  the measurements should keep track of something like this themselves
        r = calc.iterate(measurement_sweeps_per_bin)
        for p, m in six.iteritems(r):
            # NOTE: this assumes that each measurement object returns precisely a single result
            results.setdefault(p, []).append(m.get_estimate().recent_result)

    for i, m in enumerate(six.itervalues(r)):
        binlevel_data = [x for x in m.get_estimate().binlevel_data if x.nbins >= 30]
        max_error = max(x.error for x in binlevel_data)
        error_normalization = 1 / max_error if max_error else 1
        if binlevel_data[0].error:
            tau = .5 * ((max_error / binlevel_data[0].error) ** 2 - 1)
        else:
            tau = 0
        logger.info("Measurement %d relative bin errors (tau=%.2f):\t(%s)", i, tau,
                    ', '.join("{:.3f}".format(d.error * error_normalization)
                              for d in binlevel_data))

    ri = calc.simulations[0].run_information
    logger.info('compiled with "%s" (eigen %s, boost %s)', ri.compiler, ri.eigen_version, ri.boost_version)
    logger.info("%d digits of precision, exponent range [%d, %d]", ri.precision.digits,
                ri.precision.min_exponent, ri.precision.max_exponent)

    for sim in calc.simulations:
        logger.info("%s had %.2f%% of steps accepted (with %.2f%% fully rejected)",
                    sim.walk_plan.__class__.__name__,
                    (100.0 * sim.steps_accepted / sim.steps_completed),
                    (100.0 * sim.steps_fully_rejected / sim.steps_completed))

    return {p: numpy.array(r) for p, r in six.iteritems(results)}
