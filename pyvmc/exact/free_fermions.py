import numpy

from pyvmc.core import Subsystem, Orbitals, MomentaOrbitals
from pyvmc.core.orbitals import allowed_momentum

#def G_infinite(i, filling=.5):
#    if i == 0:
#        return 0.5
#    i = numpy.abs(i)
#    return numpy.sin(pi * filling * i) / pi / i

def G(lattice, orbitals):
    assert lattice.basis_indices == 1
    assert isinstance(orbitals, MomentaOrbitals)

    N_sites = len(lattice)
    green = {}
    for site in lattice:
        # FIXME: check logic of dividing by N_sites if in 2D
        green[site.bs] = sum(numpy.exp(2j * numpy.pi * numpy.dot(site.bs, allowed_momentum(q, lattice, orbitals.boundary_conditions)))
                             for q in orbitals.momentum_sites) / N_sites
    def GN(bravais_site):
        return green[bravais_site]

    return GN

def exact_renyi(subsystem, orbitals, boundary_conditions, alpha):
    assert isinstance(subsystem, Subsystem)
    subsystem_size = len(subsystem)
    lattice = subsystem.lattice
    orbitals = Orbitals.from_description(orbitals, lattice)

    corrmat = numpy.zeros(shape=(subsystem_size, subsystem_size), dtype=complex)
    GN = G(lattice, orbitals)
    assert lattice.basis_indices == 1
    for i_, i in enumerate(subsystem):
        for j_, j in enumerate(subsystem):
            site, phase = lattice.enforce_boundary(tuple(numpy.subtract(i.bs, j.bs)), boundary_conditions)
            # fixme: check phase multiply logic
            corrmat[i_][j_] = phase * GN(site)
    sign, logprefactor = numpy.linalg.slogdet(numpy.identity(subsystem_size) - corrmat)
    eigenvalues = [a / (1 - a) for a in numpy.linalg.eigvalsh(corrmat)]
    return -alpha * logprefactor - sum(numpy.log(1 + a.real ** alpha) for a in eigenvalues)
