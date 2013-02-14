#!/usr/bin/env python

import logging

import numpy

from pyvmc.core import Lattice, Bands, FreeFermionWavefunction, SimpleSubsystem, periodic
from pyvmc.core.universe import SimulationUniverse
from pyvmc.library.renyi import RenyiEntropyMeasurementPlan # RenyiLengthScalingMeasurementPlan
from pyvmc.exact import free_fermions

logger = logging.getLogger(__name__)

def test_1d_free_fermion_renyi(tolerance=0.02):
    N = 10
    F = 5

    lattice = Lattice([N])
    orbitals = Bands([F], [periodic])

    wf = FreeFermionWavefunction(lattice=lattice, orbitals=[orbitals])

    exact_values = [free_fermions.exact_renyi(SimpleSubsystem([i], lattice), orbitals, [periodic], 2)
                    for i in range(1, N // 2 + 1)]

    plans = [RenyiEntropyMeasurementPlan(wf, SimpleSubsystem([i], lattice))
             for i in range(1, N // 2 + 1)]
    calc = SimulationUniverse(plans, 500000)
    while True:
        results = calc.iterate(200000)
        measured_values = [plan.calculate(lambda p, k=None: results[p].get_estimate(k).result)
                           for plan in plans]

        differences = [measured_value - exact_value
                       for measured_value, exact_value in zip(measured_values, exact_values)]
        logger.info("differences from expected: %s", differences)
        if all(numpy.abs(d) < tolerance for d in differences):
            break

def test_2d_free_fermion_renyi(tolerance=0.02):
    lattice = Lattice([8, 8])

    # fixme: abstract away choosing the lowest F energy levels; have this
    # function return number of degenerate points at the fermi level

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_1d_free_fermion_renyi(.001)
