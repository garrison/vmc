import logging

import numpy

from pyvmc.control.scheduler import default_scheduler
from pyvmc.tmp.goals import RenyiLengthScaling
from pyvmc.core import Lattice, Bands, FreeFermionWavefunction, SimpleSubsystem, periodic, antiperiodic
from pyvmc.exact import free_fermions

def do_test_of_plan(plan, exact, get_measured, tolerance):
    def my_cb(args, plan):
        differences = [measured_value - exact_value
                       for measured_value, exact_value in zip(get_measured(plan), exact)]
        logging.info("differences from expected: %s", differences)
        if all(numpy.abs(d) < tolerance for d in differences):
            raise StopIteration

    submit_plan(plan, my_cb)
    default_scheduler.use_cores()
    default_scheduler.run()

def submit_plan(plan, callback):

    def _cb(args, plan):
        try:
            callback(args, plan)
        except StopIteration:
            from twisted.internet import reactor
            reactor.stop()
        else:
            d = plan.advance()
            d.addCallback(_cb, plan)

    d = plan.advance()
    d.addCallback(_cb, plan)

def test_1d_free_fermion_renyi(tolerance=0.02):
    N = 10
    F = 5

    lattice = Lattice([N])
    orbitals = Bands([F], [periodic])

    free_fermion = FreeFermionWavefunction(lattice=lattice, orbitals=orbitals)

    exact = [free_fermions.exact_renyi(SimpleSubsystem([i], lattice), orbitals, [periodic], 2)
             for i in xrange(1, N // 2 + 1)]

    plan = RenyiLengthScaling(free_fermion, independent=1)
    do_test_of_plan(plan, exact, lambda plan: [p.get_renyi(0) for p in plan.renyi_plans], tolerance)

def test_2d_free_fermion_renyi(tolerance=0.02):
    lattice = Lattice([8, 8])

    # fixme: abstract away choosing the lowest F energy levels; have this
    # function return number of degenerate points at the fermi level

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_1d_free_fermion_renyi(.001)
