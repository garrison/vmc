from pyvmc.core.boundary_conditions import periodic, antiperiodic

def test_boundary_condition_definitions():
    assert periodic == 1
    assert antiperiodic * 2 == 1
