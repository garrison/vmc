from pyvmc.core.boundary_conditions import periodic, antiperiodic, open_bc, periodic_bc, antiperiodic_bc

def test_boundary_condition_definitions():
    assert periodic_bc.p == periodic.p == 1
    assert antiperiodic_bc.p * 2 == antiperiodic.p * 2 == 1
    assert open_bc.p == 0

    assert periodic_bc.phase == 1
    assert antiperiodic_bc.phase == -1
    assert open_bc.phase == 0

def test_boundary_condition_as_string():
    assert str(periodic_bc) == 'periodic'
    assert str(antiperiodic_bc) == 'antiperiodic'
    assert str(open_bc) == 'open'

def test_boundary_condition_repr():
    assert repr(periodic_bc) == 'periodic_bc'
    assert repr(antiperiodic_bc) == 'antiperiodic_bc'
    assert repr(open_bc) == 'open_bc'
