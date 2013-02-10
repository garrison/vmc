#!/usr/bin/env python

import logging

import numpy

from pyvmc.core import Lattice, Bands, FreeFermionWavefunction, SimpleSubsystem, periodic
from pyvmc.exact import free_fermions

logger = logging.getLogger(__name__)

def test_1d_free_fermion_renyi(tolerance=0.02):
    N = 10
    F = 5

    lattice = Lattice([N])
    orbitals = Bands([F], [periodic])

    free_fermion = FreeFermionWavefunction(lattice=lattice, orbitals=[orbitals])

    exact = [free_fermions.exact_renyi(SimpleSubsystem([i], lattice), orbitals, [periodic], 2)
             for i in xrange(1, N // 2 + 1)]

    #plan = RenyiLengthScaling(free_fermion, independent=1)

    #differences = [measured_value - exact_value
    #               for measured_value, exact_value in zip(get_measured(plan), exact)]
    #logger.info("differences from expected: %s", differences)
    #assert all(numpy.abs(d) < tolerance for d in differences)

def test_2d_free_fermion_renyi(tolerance=0.02):
    lattice = Lattice([8, 8])

    # fixme: abstract away choosing the lowest F energy levels; have this
    # function return number of degenerate points at the fermi level

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_1d_free_fermion_renyi(.001)
