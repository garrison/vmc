import collections

import numpy
from numpy import exp, pi

from pyvmc.core.lattice import Lattice
from pyvmc.core.orbitals import allowed_momentum
from pyvmc.core.boundary_conditions import valid_boundary_conditions

# fixme: it would be much smarter to use n-dimensional fft somehow

two_pi_i = complex(0, pi + pi)

def fourier_transform(values, lattice, boundary_conditions):
    """Performs a Fourier transform on the lattice.

    `values` should be a list representing the value on each site, ordered by
    the index of each Bravais site on the lattice
    """
    assert isinstance(lattice, Lattice)
    assert isinstance(values, collections.Sequence)
    if not len(values) == lattice.total_bravais_sites():
        raise RuntimeError("wrong number of values provided for lattice fourier transform")
    assert boundary_conditions is not None
    assert valid_boundary_conditions(boundary_conditions, len(lattice.dimensions))

    return [sum([exp(-two_pi_i * numpy.dot(r, allowed_momentum(momentum_site, lattice, boundary_conditions))) * v
                 for r, v in zip(lattice.iterate_bravais_sites(), values)])
            for momentum_site in lattice.iterate_bravais_sites()]

def inverse_fourier_transform(values, lattice, boundary_conditions):
    assert isinstance(lattice, Lattice)
    assert isinstance(values, collections.Sequence)
    if not len(values) == lattice.total_bravais_sites():
        raise RuntimeError("wrong number of values provided for lattice fourier transform")
    assert boundary_conditions is not None
    assert valid_boundary_conditions(boundary_conditions, len(lattice.dimensions))

    normalization = 1.0 / len(values)

    return [sum([exp(two_pi_i * numpy.dot(r, allowed_momentum(momentum_site, lattice, boundary_conditions))) * v
                 for momentum_site, v in zip(lattice.iterate_bravais_sites(), values)]) * normalization
            for r in lattice.iterate_bravais_sites()]

