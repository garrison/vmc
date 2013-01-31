from pyvmc.core.boundary_conditions import periodic, antiperiodic, open_bc, periodic_bc, antiperiodic_bc

def test_boundary_condition_definitions():
    assert periodic_bc == periodic == 1
    assert antiperiodic_bc * 2 == antiperiodic * 2 == 1
    assert open_bc == 0
