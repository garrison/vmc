from pyvmc.core.boundary_conditions import periodic, antiperiodic, open_bc, periodic_bc, antiperiodic_bc, boundary_condition_to_string

def test_boundary_condition_definitions():
    assert periodic_bc == periodic == 1
    assert antiperiodic_bc * 2 == antiperiodic * 2 == 1
    assert open_bc == 0

def test_boundary_condition_to_string():
    assert boundary_condition_to_string(periodic) == 'periodic'
    assert boundary_condition_to_string(antiperiodic) == 'antiperiodic'
    assert boundary_condition_to_string(open_bc) == 'open'
