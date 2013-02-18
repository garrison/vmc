from pyvmc.core.lattice import LatticeSite, Lattice, HypercubicLattice, HexagonalLattice
from pyvmc.core.subsystem import Subsystem, SimpleSubsystem
from pyvmc.core.boundary_conditions import periodic, antiperiodic, periodic_bc, antiperiodic_bc, open_bc, valid_boundary_conditions, valid_closed_boundary_conditions
from pyvmc.core.orbitals import Orbitals, OrbitalsDescription, MomentaOrbitals, Bands
from pyvmc.core.wavefunction import Wavefunction, FreeFermionWavefunction

def get_vmc_version():
    import os
    from subprocess import check_output, CalledProcessError

    original_directory = os.getcwd()
    try:
        os.chdir(os.path.dirname(__file__))
        try:
            return check_output(["git", "describe", "--always"])
        except CalledProcessError:
            return None
    finally:
        os.chdir(original_directory)

import six
if six.PY3:
    import numpy
    from collections import Hashable
    from numbers import Integral

    # FIXME: track down why we are passing around these int64's anyway
    if not isinstance(numpy.int64, Integral):
        Integral.register(numpy.int64)
    if not isinstance(numpy.int64, Hashable):
        Hashable.register(numpy.int64)

    # FIXME: numpy.ndarray is Hashable in python2 but not python3.  Honestly, I
    # don't think it should be in either...
    Hashable.register(numpy.ndarray)
