import six

from cython.operator cimport dereference as deref
from libcpp.list cimport list as stdlist
from pyvmc.includes.libcpp.memory cimport unique_ptr, shared_ptr

import logging
import collections

from pyvmc.core.rng cimport RandomNumberGenerator, CppRandomNumberGenerator
from pyvmc.core.lattice cimport Lattice, CppLattice
from pyvmc.core.walk cimport Walk, CppWalk
from pyvmc.core.walk import WalkPlan
from pyvmc.core.measurement cimport BaseMeasurement, CppBaseMeasurement
from pyvmc.core.measurement import BasicMeasurementPlan
from pyvmc.core.estimate import Estimate
from pyvmc.utils.resource_logging import log_rusage
from pyvmc.utils import custom_json as json
from pyvmc.core import get_vmc_version

cdef extern from "MetropolisSimulation.hpp":
    cdef cppclass CppMetropolisSimulation "MetropolisSimulation<probability_t>":
        CppMetropolisSimulation(unique_ptr[CppWalk], stdlist[shared_ptr[CppBaseMeasurement]], unsigned int, unique_ptr[CppRandomNumberGenerator]&) nogil except +
        void iterate(unsigned int) nogil except +
        void check_for_numerical_error() nogil except +
        unsigned int steps_completed()
        unsigned int steps_accepted()
        unsigned int steps_fully_rejected()

logger = logging.getLogger(__name__)

cdef class MetropolisSimulation(object):
    cdef unique_ptr[CppMetropolisSimulation] autoptr

    cdef object _walk_plan
    cdef object _measurement_dict
    cdef unsigned int _equilibrium_steps

    cdef str _rng_name
    cdef unsigned long _rng_seed

    cdef bint _unrecoverable_error_has_occurred

    def __init__(self, walk_plan not None, Lattice lattice not None, measurement_plans, unsigned int equilibrium_steps, RandomNumberGenerator rng not None):
        """keep in mind that the rng passed can no longer be used for other things afterwards"""
        assert rng.is_good()
        self._rng_name = rng.name
        self._rng_seed = rng.seed

        self._equilibrium_steps = equilibrium_steps

        assert isinstance(walk_plan, WalkPlan)
        self._walk_plan = walk_plan
        cdef Walk walk = walk_plan.create_walk(rng)
        assert walk.autoptr.get() is not NULL

        assert isinstance(measurement_plans, collections.Sequence)
        assert all(isinstance(mp, BasicMeasurementPlan) for mp in measurement_plans)
        assert all(mp.walk_plan == walk_plan for mp in measurement_plans)
        self._measurement_dict = {mp: mp.to_measurement() for mp in measurement_plans}
        assert len(self._measurement_dict) == len(measurement_plans)  # otherwise some measurement plans are redundant

        cdef stdlist[shared_ptr[CppBaseMeasurement]] measurement_list
        cdef BaseMeasurement measurement_
        for measurement in self._measurement_dict.values():
            measurement_ = measurement
            if not measurement_.sharedptr.get().is_valid_walk(deref(walk.autoptr)):
                raise ValueError("invalid walk/measurement/wavefunction combination")
            measurement_list.push_back(measurement_.sharedptr)

        with log_rusage(logger, "Equilibrated walk using {} steps.".format(equilibrium_steps)):
            with nogil:
                self.autoptr.reset(new CppMetropolisSimulation(walk.autoptr, measurement_list, equilibrium_steps, rng.autoptr))
                # this is optional here, but we might as well before taking measurements
                self.autoptr.get().check_for_numerical_error()

        self._unrecoverable_error_has_occurred = False

        logger.info("Now prepared to consider %d different measurement(s).", len(measurement_plans))

    def iterate(self, unsigned int sweeps, check_for_numerical_error=True):
        if self._unrecoverable_error_has_occurred:
            raise RuntimeError("Unrecoverable error occurred previously during this simulation.")

        cdef bint check_for_numerical_error_ = bool(check_for_numerical_error)
        try:
            with log_rusage(logger, "Performed {} sweeps on walk.".format(sweeps)):
                with nogil:
                    self.autoptr.get().iterate(sweeps)
                    if check_for_numerical_error_:
                        self.autoptr.get().check_for_numerical_error()
        except Exception:
            self._unrecoverable_error_has_occurred = True
            raise

    def to_hdf5(self, *args, **kwargs):
        if self._unrecoverable_error_has_occurred:
            raise RuntimeError("Unrecoverable error occurred previously during this simulation.")
        return _save_simulation_to_hdf5(self, *args, **kwargs)

    # XXX: Cython causes a segfault if I declare this a staticmethod
    # (see http://trac.cython.org/cython_trac/ticket/804)
    @classmethod
    def from_hdf5(cls, *args, **kwargs):
        return RestoredSimulation(*args, **kwargs)

    property walk_plan:
        def __get__(self):
            return self._walk_plan

    property measurement_dict:
        def __get__(self):
            if self._unrecoverable_error_has_occurred:
                raise RuntimeError("Unrecoverable error occurred previously during this simulation.")
            return self._measurement_dict

    property equilibrium_steps:
        def __get__(self):
            return self._equilibrium_steps

    property steps_completed:
        def __get__(self):
            return self.autoptr.get().steps_completed()

    property steps_accepted:
        def __get__(self):
            return self.autoptr.get().steps_accepted()

    property steps_fully_rejected:
        def __get__(self):
            return self.autoptr.get().steps_fully_rejected()

    property run_information:
        def __get__(self):
            return RunInformation()

    property rng_name:
        def __get__(self):
            return self._rng_name

    property rng_seed:
        def __get__(self):
            return self._rng_seed

cdef extern from "RunInformation.hpp" namespace "RunInformation":
    const char *compiler
    const char *boost_version
    const char *eigen_version

    int precision_digits "RunInformation::Precision::digits"
    int precision_min_exponent "RunInformation::Precision::min_exponent"
    int precision_max_exponent "RunInformation::Precision::max_exponent"

cdef class RunInformation(object):
    property compiler:
        def __get__(self):
            return compiler

    property boost_version:
        def __get__(self):
            return boost_version

    property eigen_version:
        def __get__(self):
            return eigen_version

    property precision:
        def __get__(self):
            return PrecisionInformation()

cdef class PrecisionInformation(object):
    property digits:
        def __get__(self):
            return precision_digits

    property min_exponent:
        def __get__(self):
            return precision_min_exponent

    property max_exponent:
        def __get__(self):
            return precision_max_exponent

def _save_simulation_to_hdf5(sim, h5group, wf, i):
    walk_group = h5group.create_group("{0}_{1}".format(i, sim.walk_plan.__class__.__name__))

    # save the current date/time
    from datetime import datetime
    walk_group.attrs["datetime"] = str(datetime.now())

    # save information about the compilation/version
    walk_group.attrs["vmc_version"] = str(get_vmc_version())
    ri = sim.run_information
    walk_group.attrs["eigen_version"] = ri.eigen_version
    walk_group.attrs["boost_version"] = ri.boost_version
    walk_group.attrs["compiler"] = ri.compiler
    walk_group.attrs["precision_digits"] = ri.precision.digits
    walk_group.attrs["precision_min_exponent"] = ri.precision.min_exponent
    walk_group.attrs["precision_max_exponent"] = ri.precision.max_exponent

    # save rng seed and stats from the walk
    walk_group.attrs["rng_name"] = sim.rng_name
    walk_group.attrs["rng_seed"] = str(sim.rng_seed)
    walk_group.attrs["steps_accepted"] = sim.steps_accepted
    walk_group.attrs["steps_completed"] = sim.steps_completed
    walk_group.attrs["steps_fully_rejected"] = sim.steps_fully_rejected
    walk_group.attrs["equilibrium_steps"] = sim.equilibrium_steps

    # save a description of the walk that was performed
    walk_group.attrs["walkplan_json"] = json.dumps(sim.walk_plan.to_json())

    # save each measurement to an hdf5 subgroup
    for j, measurement_plan in enumerate(sim.measurement_dict):
        _save_measurement_to_hdf5(measurement_plan, walk_group, wf, sim, j)

class RestoredSimulation(object):
    def __init__(self, walk_group, wf):
        # fixme: do we want to load the date/time?
        # fixme: should we be saving the utime, stime, walltime, etc?
        # fixme: load what we now call the "run information"

        self.rng_name = walk_group.attrs["rng_name"]
        self.rng_seed = int(walk_group.attrs["rng_seed"])
        self.steps_accepted = walk_group.attrs["steps_accepted"]
        self.steps_completed = walk_group.attrs["steps_completed"]
        self.steps_fully_rejected = walk_group.attrs["steps_fully_rejected"]
        self.equilibrium_steps = walk_group.attrs["equilibrium_steps"]

        self.walk_plan = WalkPlan.from_json(json.loads(walk_group.attrs["walkplan_json"]), wf)

        self.measurement_dict = {m._measurement_plan: m
                                 for m in (RestoredMeasurement(walk_group[k], wf)
                                           for k in walk_group if k.startswith("measurement:"))}

def _save_measurement_to_hdf5(measurement_plan, walk_group, wf, sim, j):
    meas_group = walk_group.create_group("measurement:{0}_{1}".format(j, measurement_plan.__class__.__name__))

    # save a description of the measurement
    meas_group.attrs["measurementplan_json"] = json.dumps(measurement_plan.to_json())

    # save each estimate
    for key, estimate in six.iteritems(sim.measurement_dict[measurement_plan].get_estimates()):
        if key is None:
            estimate_group = meas_group
        else:
            json_key = json.dumps(key)
            if '/' in json_key:
                raise Exception("keys currently cannot contain slashes, as they represent path delimiters in hdf5")
            estimate_group = meas_group.create_group("estimate:{}".format(json_key))

        estimate.to_hdf5(estimate_group)

class RestoredMeasurement(object):
    def __init__(self, meas_group, wf):
        self._measurement_plan = BasicMeasurementPlan.from_json(json.loads(meas_group.attrs["measurementplan_json"]), wf)

        self._estimate_dict = {json.tuplize(json.loads(key.partition(':')[2])): Estimate.from_hdf5(meas_group[key])
                               for key in meas_group if key.startswith("estimate:")}
        if "result" in meas_group:
            self._estimate_dict[None] = Estimate.from_hdf5(meas_group)

    def get_estimate(self, key=None):
        return self._estimate_dict[key]

    def get_estimates(self):
        # NOTE: the caller should not modify the returned dict!
        return self._estimate_dict
