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

class SimulationUniverse(object):
    """Given a collection of measurement plans, this class will manage simulating them.
    """

    def __init__(self, measurement_plans, equilibrium_sweeps=500000):
        # NOTE: it may one day be useful to keep a list of all
        # measurement_plans originally passed, but at the moment this does not
        # seem to be useful.

        # first create all the basic measurements that need to be performed
        assert isinstance(measurement_plans, Sequence)
        assert all(isinstance(mp, MeasurementPlan) for mp in measurement_plans)
        basic_measurement_plans = frozenset(chain.from_iterable(mp.get_measurement_plans() for mp in measurement_plans))
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
            sim.iterate(sweeps)

    def get_overall_measurement_dict(self):
        return dict(chain.from_iterable(six.iteritems(sim.measurement_dict) for sim in self.simulations))

def do_calculate_plans(plans, equilibrium_sweeps=500000, bins=100, measurement_sweeps_per_bin=10000):
    # this is deprecated because it only works with measurements that have one
    # estimate.  (for instance, it does not work with
    # SubsystemOccupationProbabilityMeasurement.)
    #
    # NOTE: the bins parameter is taken into account for the number of total
    # measurement sweeps, but is currently ignored

    import warnings
    warnings.warn("do_calculate_plans() is deprecated. use SimulationUniverse directly.", DeprecationWarning, stacklevel=2)

    measurement_sweeps = bins * measurement_sweeps_per_bin

    # prepare and equilibrate simulations
    calc = SimulationUniverse(plans, equilibrium_sweeps)

    # perform simulations
    calc.iterate(measurement_sweeps)

    # log the error at each binning level
    for i, m in enumerate(six.itervalues(calc.get_overall_measurement_dict())):
        # NOTE: the following line assumes each measurement has just one estimate
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

    # log the run information
    ri = calc.simulations[0].run_information
    logger.info('compiled with "%s" (eigen %s, boost %s)', ri.compiler, ri.eigen_version, ri.boost_version)
    logger.info("%d digits of precision, exponent range [%d, %d]", ri.precision.digits,
                ri.precision.min_exponent, ri.precision.max_exponent)

    # log the acceptance rate of each simulation
    for sim in calc.simulations:
        logger.info("%s had %.2f%% of steps accepted (with %.2f%% fully rejected)",
                    sim.walk_plan.__class__.__name__,
                    (100.0 * sim.steps_accepted / sim.steps_completed),
                    (100.0 * sim.steps_fully_rejected / sim.steps_completed))

    return {mp: numpy.array(m.get_estimate().block_averages) for mp, m in six.iteritems(calc.get_overall_measurement_dict())}
