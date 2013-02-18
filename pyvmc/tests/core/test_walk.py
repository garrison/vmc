from pyvmc.core import periodic_bc
from pyvmc.core.orbitals import Bands
from pyvmc.core.lattice import Lattice
from pyvmc.core.wavefunction import FreeFermionWavefunction
from pyvmc.core.walk import WalkPlan, StandardWalkPlan

def test_walk_plan_json():
    lattice = Lattice([24, 2])
    wf = FreeFermionWavefunction(lattice=lattice, orbitals=[Bands([8, 4], [periodic_bc, periodic_bc])])
    walk_plan = StandardWalkPlan(wf)
    assert WalkPlan.from_json(walk_plan.to_json(), wf) == walk_plan
    assert StandardWalkPlan.from_json(walk_plan.to_json(), wf) == walk_plan
