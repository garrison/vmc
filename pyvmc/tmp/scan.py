from itertools import chain

import numpy

from pyvmc.utils import custom_json as json
from pyvmc.core.simulation import MetropolisSimulation
from pyvmc.core.rng import RandomNumberGenerator

def _create_universe_set(plans):
    return set(chain.from_iterable(p.get_measurement_plans() for p in plans))

def do_calculate_plans(plans):
    # first get all the measurements that need to be performed
    universe = {p: p.to_measurement() for p in _create_universe_set(plans)}

    # now organize them by each walk which must be performed
    by_walk = {}
    for p, m in universe.iteritems():
        by_walk.setdefault(p.walk, []).append(m)

    # prepare and equilibriate simulations
    sims = []
    for walk, measurements in by_walk.iteritems():
        sims.append(MetropolisSimulation(walk.create_walk(RandomNumberGenerator()), walk.wavefunction.lattice, measurements, 500000))

    # perform simulations
    universe_results = {p: [] for p in universe}
    for i in xrange(100):
        for sim in sims:
            sim.iterate(10000)
        for p, m in universe.iteritems():
            # fixme: can we guarantee that p.to_json() doesn't have any forward slashes?
            # fixme: we should really save each measurement in a subgroup for that walk, and in there store information from the walk's completion... INCLUDING rusage, etc

            universe_results[p].append(m.get_result())
            # FIXME: do a reset!

    return universe_results

def calculate_plans(plan_dict, h5group):
    universe_results = do_calculate_plans(plan_dict.values())

    # now save the results
    for p, results in universe_results.iteritems():
        dataset = h5group.create_dataset(json.dumps(p.to_json()), data=numpy.array(results))
        dataset.attrs["NOTE"] = "measurement reset not yet implemented"
    h5group.file.flush()

class ResultReturner(object):
    def __init__(self, result):
        self.result = result

    def get_result(self):
        return self.result

def load_results(plan_dict, h5group):
    universe_set = _create_universe_set(plan_dict.values())
    universe = {plan: ResultReturner(h5group[json.dumps(plan.to_json())][-1]) for plan in universe_set}
    results = {k: plan.get_result(universe) for k, plan in plan_dict.iteritems()}

    return results
