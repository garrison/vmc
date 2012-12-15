from pyvmc.core.lattice import Lattice
from pyvmc.core.boundary_conditions import periodic, antiperiodic
from pyvmc.core.fourier import fourier_transform, inverse_fourier_transform

def test_fourier():
    lattice = Lattice([2, 2])
    boundary_conditions = (periodic, antiperiodic)

    x = [0.3, 0.76, 0.1j, 0.83]
    y = fourier_transform(x, lattice, boundary_conditions)

    # test inverse fourier transform
    z = inverse_fourier_transform(y, lattice, boundary_conditions)
    err = sum(abs(a - b) for a, b in zip(x, z))
    assert err < 1e-10
