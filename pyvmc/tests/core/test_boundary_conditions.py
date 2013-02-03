from pyvmc.core.boundary_conditions import periodic, antiperiodic, open_bc, periodic_bc, antiperiodic_bc, boundary_condition_to_string

def test_boundary_condition_definitions():
    assert periodic_bc.p == periodic.p == 1
    assert antiperiodic_bc.p * 2 == antiperiodic.p * 2 == 1
    assert open_bc.p == 0

    assert periodic_bc.phase == 1
    assert antiperiodic_bc.phase == -1
    assert open_bc.phase == 0

def test_boundary_condition_to_string():
    assert boundary_condition_to_string(periodic) == 'periodic'
    assert boundary_condition_to_string(antiperiodic) == 'antiperiodic'
    assert boundary_condition_to_string(open_bc) == 'open'
