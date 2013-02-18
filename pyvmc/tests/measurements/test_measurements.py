import pytest

from pyvmc.core import periodic_bc
from pyvmc.core.orbitals import Bands
from pyvmc.core.lattice import Lattice
from pyvmc.core.subsystem import SimpleSubsystem
from pyvmc.core.wavefunction import FreeFermionWavefunction
from pyvmc.core.measurement import BasicMeasurementPlan
from pyvmc.measurements import SubsystemOccupationProbabilityMeasurementPlan

def test_measurement_plan_json():
    lattice = Lattice([24, 2])
    wf = FreeFermionWavefunction(lattice=lattice, orbitals=[Bands([8, 4], [periodic_bc, periodic_bc])])
    subsystem = SimpleSubsystem([8, 2], lattice)

    sonpm_plan = SubsystemOccupationProbabilityMeasurementPlan(wf, subsystem, steps_per_measurement=500)
    assert BasicMeasurementPlan.from_json(sonpm_plan.to_json(), wf) == sonpm_plan
    assert SubsystemOccupationProbabilityMeasurementPlan.from_json(sonpm_plan.to_json(), wf) == sonpm_plan
