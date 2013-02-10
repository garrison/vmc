from __future__ import division

# fixme: remove need for copy, deepcopy, and hashable_json

from pyvmc.utils import custom_json as json
import random
from copy import deepcopy, copy
import logging

from twisted.internet import defer, task, reactor

## fixme: or we pass a scheduler in advance()
#from pyvmc.control.scheduler import default_scheduler

from pyvmc.core.simulation import MetropolisSimulation
from pyvmc.core.rng import RandomNumberGenerator

logger = logging.getLogger(__name__)

def hashable_json(obj):
    if isinstance(obj, dict):
        return frozenset((k, hashable_json(v)) for k, v in obj.iteritems())
    elif isinstance(obj, list):
        return tuple(hashable_json(a) for a in obj)
    else:
        return obj

class TmpMeasurementPlan(object):
    def __init__(self, measurement_plan, parse_result):
        self.lattice = measurement_plan.walk.wavefunction.lattice
        self.walk = {}
        self.walk["system"] = deepcopy(measurement_plan.walk.wavefunction.to_json())
        self.walk["simulation"] = deepcopy(measurement_plan.walk.to_json())
        self.measurement_plan = measurement_plan
        self.hashable = hashable_json([self.walk, measurement_plan])
        self.hashable_walk = hashable_json(self.walk)
        self.parse_result = parse_result  # fixme (?)

class Walk(object):
    equilibrium = False
    sim = None
    measurement_steps_completed = 0

    def __init__(self, walk_json, ww):
        self.walk_json = deepcopy(walk_json)
        self.ww = ww
        self.equilibrium_steps = 500000
        self.measurements_in_progress = []
        self.measurements_pending = []
        self.waiting_on_walk = False

    def _advancement_completed(self, args):
        logger.info("monte carlo stats: %s", {
            "steps-completed": self.sim.steps_completed,
            "steps-accepted": self.sim.steps_accepted,
            "steps-fully-rejected": self.sim.steps_fully_rejected,
        })
        for measurement, deferred in self.measurements_in_progress:
            estimates = measurement.measurement.get_estimates()
            try:
                r = estimates[None].result
            except KeyError:
                # it must have returned a full dictionary of results
                r = [(k, v.result) for k, v in estimates.items()]
            deferred.callback((r,))
        del self.measurements_in_progress[:]
        if self.measurements_pending:
            reactor.callLater(0, self._advance_pending_measurements)
            assert not self.waiting_on_walk
            self.waiting_on_walk = True

    def _advancement_failed(self, args):
        for measurement, deferred in self.measurements_in_progress:
            deferred.errback()
        del self.measurements_in_progress[:]

    def _advance_pending_measurements(self):
        self.waiting_on_walk = False
        self.measurements_in_progress = self.measurements_pending
        self.measurements_pending = []
        try:
            if self.sim is None:
                # XXX FIXME: if measurements_in_progress changes between
                # advances, this WILL NOT notice and will continue with the
                # old measurements.
                rng = RandomNumberGenerator()
                self.sim = MetropolisSimulation(self.ww.create_walk(rng),
                                                self.measurements_in_progress[0][0].measurement_plan.lattice,
                                                [m[0].measurement for m in self.measurements_in_progress],
                                                self.equilibrium_steps,
                                                rng)
            # the following will always result in the number of steps completed
            # being a power of two
            steps = self.measurement_steps_completed or 512
            self.sim.iterate(steps)
            self.measurement_steps_completed += steps
        except Exception as e:
            self._advancement_failed(None)
        else:
            self._advancement_completed(None)

    def advance_measurement(self, measurement):
        assert (measurement not in self.measurements_pending)
        deferred = defer.Deferred()
        self.measurements_pending.append((measurement, deferred))
        if not self.measurements_in_progress and not self.waiting_on_walk:
            reactor.callLater(0, self._advance_pending_measurements)
        self.waiting_on_walk = True
        return deferred

class WalkSet(object):
    def __init__(self, walk_json, ww, n_independent):
        # walk_json has already been deepcopy'd, so we can just store a
        # reference to it and it shouldn't change.
        self.walk_json = walk_json
        self.ww = ww
        self.ws = [Walk(walk_json, ww) for i in xrange(n_independent)]

    def extend_by(self, n):
        self.ms.extend(Walk(self.walk_json, self.ww) for i in xrange(n))

    def extend_to(self, n):
        if n > len(self.ws):
            self.extend_by(n - len(self.ws))

    def __len__(self):
        return len(self.ws)

    def __iter__(self):
        return iter(self.ws)

    def __getitem__(self, index):
        return self.ws[index]

class Measurement(object):
    def __init__(self, measurement_plan):
        self.measurement_plan = measurement_plan
        self.measurement = measurement_plan.measurement_plan.to_measurement()
        self.result = None
        self.has_been_advanced = False

    def _record_result(self, args):
        self._add_result(*args)
        return None

    def advance(self, walk):
        d = walk.advance_measurement(self)
        d.addCallback(self._record_result)
        return d

    def _add_result(self, result_json):
        result = self.measurement_plan.parse_result(result_json)
        self.result = copy(result)

        self.has_been_advanced = True

    def get_aggregate_result(self):
        if not self.has_been_advanced:
            # or we could just return None here...
            raise Exception("no measurements have been added.  cannot take aggregate result.")

        return copy(self.result)

class MeasurementSet(object):
    """A set of identical, actual measurements taken directly from independent monte carlo runs"""

    def __init__(self, measurement_plan, walk_set, n_independent):
        self.measurement_plan = measurement_plan
        self.walk_set = walk_set
        self.ms = [Measurement(measurement_plan) for i in xrange(n_independent)]
        self.advancement_in_progress_deferreds = []  # blank if not advancement is in progress

    def extend_by(self, n):
        self.ms.extend(Measurement(self.measurement_plan) for i in xrange(n))
        self.walk_set.extend_to(len(self.ms))

    def extend_to(self, n):
        if n > len(self.ms):
            self.extend_by(n - len(self.ms))

    def __len__(self):
        return len(self.ms)

    def __iter__(self):
        return iter(self.ms)

    def __getitem__(self, index):
        return self.ms[index]

    def _advancement_completed(self, args):
        for deferred in self.advancement_in_progress_deferreds:
            deferred.callback(None)
        del self.advancement_in_progress_deferreds[:]

    def _advancement_failed(self, failure):
        for deferred in self.advancement_in_progress_deferreds:
            deferred.errback(failure)
        del self.advancement_in_progress_deferreds[:]

    def _advance(self):
        deferreds = [m.advance(self.walk_set[i]) for i, m in enumerate(self.ms)]
        return defer.gatherResults(deferreds, consumeErrors=True)

    def advance(self):
        deferred = defer.Deferred()
        if not self.advancement_in_progress_deferreds:
            self._advance().addCallbacks(self._advancement_completed,
                                         self._advancement_failed)
        self.advancement_in_progress_deferreds.append(deferred)
        return deferred

class Universe(object):
    def __init__(self):
        # fixme: i bet these could even be weakref dictionaries
        self.walksets = {}
        self.measurementsets = {}

    def _get_walk_set(self, measurement_plan, n_independent):
        try:
            walk_set = self.walksets[measurement_plan.hashable_walk]
        except KeyError:
            walk_set = WalkSet(measurement_plan.walk, measurement_plan.measurement_plan.walk, n_independent)
            self.walksets[measurement_plan.hashable_walk] = walk_set
        else:
            walk_set.extend_to(n_independent)
        return walk_set

    def get_measurement_set(self, measurement_plan, n_independent):
        try:
            measurement_set = self.measurementsets[measurement_plan.hashable]
        except KeyError:
            walk_set = self._get_walk_set(measurement_plan, n_independent)
            measurement_set = MeasurementSet(measurement_plan, walk_set, n_independent)
            self.measurementsets[measurement_plan.hashable] = measurement_set
        else:
            measurement_set.extend_to(n_independent)
        return measurement_set

global_universe = Universe()

def get_default_universe():
    return global_universe
