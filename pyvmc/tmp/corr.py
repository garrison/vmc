import sys
import itertools

import numpy
import h5py

from pyvmc.core.measurement import BasicOperatorMeasurementPlan
from pyvmc.core.operator import SiteHop, BasicOperator
from pyvmc.core import LatticeSite, periodic, antiperiodic
from pyvmc.core.simulation import MetropolisSimulation
from pyvmc.core.rng import RandomNumberGenerator

def calculate_correlators(wf, filename):
    i = LatticeSite((0, 0))
    plans = [
        (
            # green
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, j, 0)], True, (periodic, periodic))),
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, j, 1)], True, (periodic, periodic))),
            # density-density terms
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, i, 0), SiteHop(j, j, 1)], True, (periodic, periodic))),
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, i, 1), SiteHop(j, j, 0)], True, (periodic, periodic))),
            # terms in density-density and spin-spin
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, i, 0)] + ([SiteHop(j, j, 0)] if j != i else []), True, (periodic, periodic))),
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, i, 1)] + ([SiteHop(j, j, 1)] if j != i else []), True, (periodic, periodic))),
            # spin-spin terms
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, j, 0), SiteHop(j, i, 1)], True, (periodic, periodic))),
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, j, 1), SiteHop(j, i, 0)], True, (periodic, periodic))),
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, i, 0), SiteHop(j, j, 1)], True, (periodic, periodic))),
            BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(i, i, 1), SiteHop(j, j, 0)], True, (periodic, periodic))),
        )
        for j in wf.lattice
    ]
    measurements = [[plan.to_measurement() for plan in pp] for pp in plans]

    more_plans = [
        # ring-exchange terms
        BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 1)], True, (periodic, periodic))),
        BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((1, 0)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((0, 1)), 0)], True, (periodic, periodic))),
        BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 0), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 1)], True, (periodic, periodic))),
        BasicOperatorMeasurementPlan(wf, BasicOperator([SiteHop(LatticeSite((0, 0)), LatticeSite((0, 1)), 1), SiteHop(LatticeSite((1, 1)), LatticeSite((1, 0)), 0)], True, (periodic, periodic))),
    ]
    more_measurements = [plan.to_measurement() for plan in more_plans]

    all_plans = list(itertools.chain.from_iterable(plans)) + more_plans
    all_measurements = list(itertools.chain.from_iterable(measurements)) + more_measurements

    walk = plans[0][0].walk
    sim = MetropolisSimulation(walk.create_walk(RandomNumberGenerator()), wf.lattice, all_measurements, 500000)

    for zzz in xrange(500):
        sim.iterate(10000)

        green0 = numpy.ndarray(shape=wf.lattice.dimensions, dtype=complex)
        green1 = numpy.ndarray(shape=wf.lattice.dimensions, dtype=complex)
        dd = numpy.ndarray(shape=wf.lattice.dimensions, dtype=complex)
        ss = numpy.ndarray(shape=wf.lattice.dimensions, dtype=complex)
        for site, measurement in zip(wf.lattice, measurements):
            green0[site.bs] = measurement[0].get_result() / len(wf.lattice)
            green1[site.bs] = measurement[1].get_result() / len(wf.lattice)
            dd[site.bs] = sum([measurement[i].get_result() for i in xrange(2, 6)]) / len(wf.lattice) - (wf.rho ** 2)
            ss[site.bs] = (-.5 * sum([measurement[i].get_result() for i in xrange(6, 8)])
                           + .25 * sum([measurement[i].get_result() for i in xrange(4, 6)])
                           - .25 * sum([measurement[i].get_result() for i in xrange(8, 10)])
            ) / len(wf.lattice)
        ring = sum([more_measurements[i].get_result() for i in xrange(4)]).real / len(wf.lattice)  # one half of four terms plus their hermitian conjugates, average per site
        ss[(0, 0)] = .75 * wf.rho  # same-site commutation relations lead to a different result
        green0_fourier = numpy.fft.fftn(green0)
        green1_fourier = numpy.fft.fftn(green1)
        dd_fourier = numpy.fft.fftn(dd)
        ss_fourier = numpy.fft.fftn(ss)

        # save data
        with h5py.File(filename, "w") as f:
            f.attrs["wf"] = repr(wf)
            # or require_dataset(exact=True)
            f.create_dataset("Green_0", data=green0)
            f.create_dataset("GreenFourier_0", data=green0_fourier)
            f.create_dataset("Green_1", data=green1)
            f.create_dataset("GreenFourier_1", data=green1_fourier)
            f.create_dataset("DensityDensity", data=dd)
            f.create_dataset("DensityDensityFourier", data=dd_fourier)
            f.create_dataset("SpinSpin", data=ss)
            f.create_dataset("SpinSpinFourier", data=ss_fourier)
            f.create_dataset("Ring", data=ring)
            #f.flush()

            g = f.require_group("OperatorMeasurement")
            om_data = {}
            for plan, measurement in zip(all_plans, all_measurements):
                om_data[repr((plan.hops, plan.sum, plan.boundary_conditions))] = measurement.get_result()
            for k, v in om_data.iteritems():
                g.create_dataset(k, data=v)

        sys.stderr.write(".")
        sys.stderr.flush()
