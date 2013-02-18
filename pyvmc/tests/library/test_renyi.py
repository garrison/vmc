import pytest

from pyvmc.core import periodic_bc
from pyvmc.core.orbitals import Bands
from pyvmc.core.lattice import Lattice
from pyvmc.core.subsystem import SimpleSubsystem
from pyvmc.core.wavefunction import FreeFermionWavefunction
from pyvmc.core.walk import WalkPlan
from pyvmc.library.renyi import RenyiModPossibleWalkPlan, RenyiSignWalkPlan

def test_walk_plan_json():
    lattice = Lattice([24, 2])
    wf = FreeFermionWavefunction(lattice=lattice, orbitals=[Bands([8, 4], [periodic_bc, periodic_bc])])
    subsystem = SimpleSubsystem([8, 2], lattice)

    renyi_modpossible_plan = RenyiModPossibleWalkPlan(wf, subsystem)
    assert WalkPlan.from_json(renyi_modpossible_plan.to_json(), wf) == renyi_modpossible_plan
    assert RenyiModPossibleWalkPlan.from_json(renyi_modpossible_plan.to_json(), wf) == renyi_modpossible_plan

    renyi_sign_plan = RenyiSignWalkPlan(wf, subsystem)
    assert WalkPlan.from_json(renyi_sign_plan.to_json(), wf) == renyi_sign_plan
    assert RenyiSignWalkPlan.from_json(renyi_sign_plan.to_json(), wf) == renyi_sign_plan

    with pytest.raises(Exception):
        RenyiSignWalkPlan.from_json(renyi_modpossible_plan.to_json(), wf)
