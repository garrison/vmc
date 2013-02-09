from itertools import chain
import logging

import numpy

from pyvmc.utils import custom_json as json
from pyvmc.core.simulation import MetropolisSimulation
from pyvmc.core.rng import RandomNumberGenerator
from pyvmc.utils import average

logger = logging.getLogger(__name__)

def _create_universe_set(plans):
    return set(chain.from_iterable(p.get_measurement_plans() for p in plans))

def do_calculate_plans(plans, equilibrium_sweeps=500000, bins=100, measurement_sweeps_per_bin=10000):
    # first get all the measurements that need to be performed
    universe = {p: p.to_measurement() for p in _create_universe_set(plans)}

    # now organize them by each walk which must be performed
    by_walk = {}
    for p, m in universe.iteritems():
        by_walk.setdefault(p.walk, []).append(m)

    # prepare and equilibriate simulations
    sims = []
    for walk, measurements in by_walk.iteritems():
        sims.append(MetropolisSimulation(walk.create_walk(RandomNumberGenerator()), walk.wavefunction.lattice, measurements,  equilibrium_sweeps))

    # perform simulations
    universe_results = {p: [] for p in universe}
    for i in xrange(bins):  # FIXME: this outer loop should be unnecessary.  the measurements should keep track of something like this themselves
        for sim in sims:
            sim.reset_measurement_estimates()
            sim.iterate(measurement_sweeps_per_bin)
        for p, m in universe.iteritems():
            # fixme: can we guarantee that p.to_json() doesn't have any forward slashes?
            # fixme: we should really save each measurement in a subgroup for that walk, and in there store information from the walk's completion... INCLUDING rusage, etc
            universe_results[p].append(m.get_estimate().recent_result)   # NOTE: this assumes that each measurement object returns precisely a single result

    for sim, walk in zip(sims, by_walk):
        logger.info("%s had %.2f%% of steps accepted (with %.2f%% fully rejected)",
                    walk.__class__.__name__,
                    (100.0 * sim.steps_accepted / sim.steps_completed),
                    (100.0 * sim.steps_fully_rejected / sim.steps_completed))

    return {p: numpy.array(results) for p, results in universe_results.iteritems()}

def calculate_plans(plan_dict, h5group):
    universe_results = do_calculate_plans(plan_dict.values())

    # now save the results
    for p, results in universe_results.iteritems():
        dataset = h5group.create_dataset(json.dumps(p.to_json()), data=results)
    h5group.file.flush()

def load_results(plan_dict, h5group):
    return {p: h5group[json.dumps(p.to_json())] for p in _create_universe_set(plan_dict.values())}

# OLD

class ResultReturner(object):
    def __init__(self, result):
        self.result = result

    def get_result(self):
        return self.result

def load_results(plan_dict, h5group):
    universe_set = _create_universe_set(plan_dict.values())
    universe = {plan: ResultReturner(average(numpy.array(h5group[json.dumps(plan.to_json())]))) for plan in universe_set}
    results = {k: plan.get_result(universe) for k, plan in plan_dict.iteritems()}

    return results
