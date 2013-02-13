#!/usr/bin/env python

from pyvmc.core import Lattice, LatticeSite, HexagonalLattice, periodic, antiperiodic
from pyvmc.core.fourier import fourier_transform, inverse_fourier_transform
from pyvmc.library.mft_utils import plot_pairing_matrix

import numpy as np
from numpy import pi
np.set_printoptions(precision=4, suppress=True, threshold=1000, linewidth=300)

import matplotlib.pyplot as plt


def triangular_lattice_dwave():
    lattice = HexagonalLattice([12, 12])
    boundary_conditions = (periodic, antiperiodic)

    #alpha = pi/100
    #t1 = np.cos(alpha)
    #delta1 = np.sin(alpha)

    t1 = 1.0
    delta1 = 5.0

    mu = -0.3

    Nsites = len(lattice)

    t_R = np.zeros(Nsites, dtype=complex)
    delta_R = np.zeros(Nsites, dtype=complex)

    n = 0
    for nn in lattice.nearest_neighbors(lattice[0], double_count=True):
        site_R, phase = lattice.enforce_boundary(nn, boundary_conditions)
        ind_R = lattice.index(site_R)
        t_R[ind_R] = t1 * phase
    #    delta_R[ind_R] = delta1 * np.exp(2*pi*1j*n/3.0) * phase  # d+id
        delta_R[ind_R] = delta1 * np.exp(2*pi*1j*n/3.0).real * phase  # d
        n = n+1

    epsilon_k = fourier_transform(list(-t_R), lattice, boundary_conditions)
    xi_k = [ep-mu for ep in epsilon_k]
    delta_k = fourier_transform(list(delta_R), lattice, boundary_conditions)

    E_k = [np.sqrt(xi*xi + np.abs(delta)*np.abs(delta)) for xi, delta in zip(xi_k, delta_k)]
    posEvals = np.sort(np.array(E_k, dtype=complex).real)

    g_k = [delta / ( xi + E ) for xi, delta, E in zip(xi_k, delta_k, E_k)]

    phi_r = inverse_fourier_transform(g_k, lattice, boundary_conditions)
    phi_r = np.array(phi_r, dtype=complex)

    print 'phi_r:'
    print phi_r
    print

    phiMat = np.zeros((Nsites, Nsites), dtype=complex)
    for ind_r in range(0, Nsites):
        for ind_rp in range(0, Nsites):
            site_r = lattice[ind_r]
            site_rp = lattice[ind_rp]
            site_r_minus_rp, phase = lattice.enforce_boundary(site_r.negative_displaced(site_rp.bs), boundary_conditions)
            ind_r_minus_rp = lattice.index(site_r_minus_rp)
            phiMat[ind_r, ind_rp] = phase * phi_r[ind_r_minus_rp]

    # Check that phiMat is symmetric:
    #assert( np.max(np.abs(phiMat - phiMat.transpose())) < 1e-12 )
    print '"Symmetricness" of phiMat = ', np.max(np.abs(phiMat - phiMat.transpose()))

    #np.savetxt('phiMat_TI.dat', phiMat.view(float))


    ##### Plot some stuff #####

    fig = plt.figure(figsize=(17,15))

    plt.plot([ xi.real for xi in xi_k], 'o-', label=r'$\xi_k$')
    plt.plot([ delta.real for delta in delta_k ], 's-', label=r'$\Delta_k$')
    plt.plot([ E.real for E in E_k ], '^-', label=r'$E_k=\sqrt{\xi_k^2 + |\Delta_k|^2}$')
    plt.plot([ (xi + E).real for xi, E in zip(xi_k, E_k) ], 'd-', label=r'$\xi_k + E_k$')

    plt.axhline(y=0, linewidth=1, color='k', linestyle='--')
    plt.xlim(0, len(xi_k))
    plt.setp(plt.gca().get_xmajorticklabels(), fontsize=25)
    plt.setp(plt.gca().get_ymajorticklabels(), fontsize=25)
    plt.xlabel('$k$', fontsize=25)
    plt.ylabel('Energies', fontsize=25)
    plt.legend(loc='best', fontsize=20)

    fig = plt.figure(figsize=(17,15))

    plt.plot([ g.real for g in g_k ], 'o-', label=r'$g_k$')
    plt.plot([ phi.real for phi in phi_r ], 's-', label=r'$\phi(r)$')

    plt.axhline(y=0, linewidth=1, color='k', linestyle='--')
    plt.xlim(0, len(phi_r))
    plt.setp(plt.gca().get_xmajorticklabels(), fontsize=25)
    plt.setp(plt.gca().get_ymajorticklabels(), fontsize=25)
    plt.xlabel('$k,\,r$', fontsize=25)
    plt.ylabel('Amplitudes', fontsize=25)
    plt.legend(loc='best', fontsize=20)

    plot_pairing_matrix(phiMat, callshow=False)

    plt.show()

if __name__ == '__main__':
    triangular_lattice_dwave()
