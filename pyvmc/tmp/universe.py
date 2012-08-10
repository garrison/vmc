from __future__ import division

# fixme: remove need for copy, deepcopy, and hashable_json

from pyvmc.utils import complex_json as json
import random
from copy import deepcopy, copy
import logging

from twisted.internet import defer, task, reactor

## fixme: or we pass a scheduler in advance()
#from pyvmc.control.scheduler import default_scheduler

from pyvmc.core.simulation import perform_simulation

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
        self.walk = {}
        self.walk["system"] = deepcopy(measurement_plan.walk.wavefunction.to_json())
        self.walk["simulation"] = deepcopy(measurement_plan.walk.to_json())
        self.measurement = deepcopy(measurement_plan.to_json())
        self.hashable = hashable_json([self.walk, measurement_plan])
        self.hashable_walk = hashable_json(self.walk)
        self.parse_result = parse_result  # fixme (?)

class Walk(object):
    equilibrium = False

    def __init__(self, walk_json):
        self.walk_json = deepcopy(walk_json)
        self.walk_json["simulation"]["equilibrium-steps"] = 500000
        self.walk_json["simulation"]["measurement-steps"] = 500000
        self.measurements_in_progress = []
        self.measurements_pending = []
        self.waiting_on_walk = False

    def _advancement_completed(self, args):
        output, err = args
        output = json.loads(output)
        logger.info("monte carlo stats: %s", output["monte-carlo-stats"])
        self.walk_json["simulation"]["equilibrium-steps"] = 0
        self.walk_json["simulation"]["initial-positions"] = output["final-positions"]
        for result_json, (measurement, deferred) in zip(output["measurements"], self.measurements_in_progress):
            deferred.callback((result_json, self.walk_json["simulation"]["measurement-steps"]))
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
        vmc_core_input = copy(self.walk_json)
        vmc_core_input["simulation"] = copy(vmc_core_input["simulation"])
        vmc_core_input["simulation"]["measurements"] = [m[0].measurement_plan.measurement
                                                        for m in self.measurements_in_progress]
        vmc_core_input["rng"] = { "seed": random.randint(0, 2 ** 32 - 1) }
        try:
            output_string = perform_simulation(json.dumps(vmc_core_input))
        except Exception as e:
            self._advancement_failed(None)
        else:
            self._advancement_completed((output_string, None))

    def advance_measurement(self, measurement):
        assert (measurement not in self.measurements_pending)
        deferred = defer.Deferred()
        self.measurements_pending.append((measurement, deferred))
        if not self.measurements_in_progress and not self.waiting_on_walk:
            reactor.callLater(0, self._advance_pending_measurements)
        self.waiting_on_walk = True
        return deferred

class WalkSet(object):
    def __init__(self, walk_json, n_independent):
        # walk_json has already been deepcopy'd, so we can just store a
        # reference to it and it shouldn't change.
        self.walk_json = walk_json
        self.ws = [Walk(walk_json) for i in xrange(n_independent)]

    def extend_by(self, n):
        self.ms.extend(Walk(self.walk_json) for i in xrange(n))

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
        self.result = None
        self.total_measurements = 0

    def _record_result(self, args):
        self._add_result(*args)
        return None

    def advance(self, walk):
        d = walk.advance_measurement(self)
        d.addCallback(self._record_result)
        return d

    def _add_result(self, result_json, n_measurements):
        result = self.measurement_plan.parse_result(result_json)
        if self.result is None:
            if isinstance(result, dict):
                selfresult = {}
                for k, v in result.iteritems():
                    selfresult[k] = v * n_measurements
                self.result = selfresult
            else:
                self.result = result * n_measurements
        else:
            if isinstance(result, dict):
                selfresult = self.result
                assert isinstance(selfresult, dict)
                keys = set(selfresult)
                keys.update(result)
                for k in keys:
                    selfresult[k] = selfresult.get(k, 0.0) + result.get(k, 0.0) * n_measurements
            else:
                self.result += result * n_measurements

        self.total_measurements += n_measurements

    def get_aggregate_result(self):
        # fixme: cache this until it changes
        selfresult = self.result
        total_measurements = self.total_measurements

        if total_measurements == 0:
            # or we could just return None here...
            raise Exception("no measurements have been added.  cannot take aggregate result.")

        if isinstance(selfresult, dict):
            rv = {}
            for k, v in selfresult.iteritems():
                # fixme: we might want to perform some transformation on the
                # keys, for instance, to turn them into floats instead of
                # indices for k vectors.  maybe a generic mechanism here?  call
                # a classmethod that by default just returns the value...
                #
                # but then again, we will probably do more on this, such as
                # error analysis based on the key, so it is probably best just
                # to keep this elsewhere.
                rv[k] = v / total_measurements
            return rv
        else:
            return selfresult / total_measurements

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
            walk_set = WalkSet(measurement_plan.walk, n_independent)
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
