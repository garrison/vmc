#!/usr/bin/env python

from pyvmc.core import Lattice, LatticeSite, HexagonalLattice, periodic, antiperiodic

import numpy
from numpy import pi
numpy.set_printoptions(precision=15, suppress=True, threshold=10000, linewidth=300)

from scipy import optimize

import logging
logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt


##### General stuff #####

def bcs_soln(tt, delta):
    assert isinstance(tt, numpy.ndarray)
    assert isinstance(delta, numpy.ndarray)
    assert tt.ndim == delta.ndim == 2
    assert tt.shape[0] == tt.shape[1] == delta.shape[0] == delta.shape[1]

    Nsites = tt.shape[0]

    H = numpy.zeros((Nsites, Nsites, 2, 2), dtype=complex)
    for ind_r in range(0, Nsites):
        for ind_rp in range(0, Nsites):
            H[ind_r, ind_rp] = -numpy.array([[tt[ind_r, ind_rp], delta[ind_r, ind_rp]],
                                             [delta[ind_r, ind_rp].conjugate(), -tt[ind_r, ind_rp].conjugate()]])

    Hflat = numpy.zeros((2*Nsites, 2*Nsites), dtype=complex)
    for i in range(0, 2*Nsites):
        for j in range(0, 2*Nsites):
            Hflat[i,j] = H[i//2, j//2, numpy.mod(i,2), numpy.mod(j,2)]

    # Check that Hflat is Hermitian:
    assert numpy.max(numpy.abs(Hflat - Hflat.conjugate().transpose())) < 1e-12

    eigsys = numpy.linalg.eigh(Hflat)

    idx = eigsys[0].argsort()
    evals = eigsys[0][idx]
    evecs = eigsys[1][:,idx]

    # Check that eigenvectors are orthonormal:
    Icheck = numpy.dot( evecs.conjugate().transpose(), evecs )
    assert numpy.max(numpy.abs(Icheck - numpy.identity(Icheck.shape[0], dtype=complex))) < 1e-12

    evalsNeg = evals[0:Nsites]
    evecsNeg = evecs[:, 0:Nsites]
    evalsPos = evals[Nsites:2*Nsites]
    evecsPos = evecs[:, Nsites:2*Nsites]

    logger.debug(' positive eigenvalues = ...')
    logger.debug(evalsPos)
#    logger.debug(' negative eigenvalues = ...')
#    logger.debug(evalsNeg[::-1])

    # Check that positive eigenvalues are positive and that eigenvalues have come in plus-minus pairs:
    assert numpy.min(evalsPos > 0)
    assert numpy.max(numpy.abs(evalsNeg[::-1] + evalsPos)) < 1e-12

    u = evecsPos[0:evecsPos.shape[0]:2, :]
    v = evecsPos[1:evecsPos.shape[0]:2, :]

    # Check that (v, -u)^* from evecsPos are indeed negative eigenvectors:
    errr = []
    for n in range(0, Nsites):
        evecTest = numpy.zeros(2*Nsites, dtype=complex)
        for i in range(0, Nsites):
            evecTest[2*i] = v[i,n].conjugate()
            evecTest[2*i+1] = -u[i,n].conjugate()
        errr.append(numpy.max(numpy.abs( numpy.dot( Hflat, evecTest ) - (-evalsPos[n] * evecTest ))))
    assert numpy.max(errr) < 1e-12

    return {'u':u, 'v':v, 'evalsPos':evalsPos}


def bcs_stats(u, v):
    assert isinstance(u, numpy.ndarray)
    assert isinstance(v, numpy.ndarray)
    assert u.ndim == v.ndim == 2
    assert u.shape[0] == u.shape[1] == v.shape[0] == v.shape[1]

    Nsites = u.shape[0]

    fdagf = numpy.zeros((Nsites), dtype=float)
    fdownfup = numpy.zeros((Nsites), dtype=complex)
    Tvec = numpy.zeros((Nsites, 3), dtype=float)
    Ttot = numpy.zeros((Nsites), dtype=float)

    for ind_r in range(0, Nsites):
        fdagf[ind_r] = 2.0 * numpy.dot( v[ind_r, :], v[ind_r, :].conjugate() ).real
        fdownfup[ind_r] = -numpy.dot( u[ind_r, :], v[ind_r, :].conjugate() )

    for ind_r in range(0, Nsites):
        Tvec[ind_r, 0] = fdownfup[ind_r].real
        Tvec[ind_r, 1] = -fdownfup[ind_r].imag
        Tvec[ind_r, 2] = 0.5 * (fdagf[ind_r] - 1.0)
        Tvec[ind_r, :] = 2.0 * Tvec[ind_r, :]  # 2.0 to be consistent with Lesik.
        Ttot[ind_r] = numpy.sqrt( numpy.dot( Tvec[ind_r, :], Tvec[ind_r, :].conjugate() ) ).real

    return {'fdagf':fdagf, 'fdownfup':fdownfup, 'Tvec':Tvec, 'Ttot':Ttot}


def bcs_stats_average(u, v):
    stats = bcs_stats(u, v)

    if numpy.max(numpy.std(stats['Tvec'], 0)) > 1e-10:
        logger.warning(' This BCS state is not translationally invariant!  standard deviation of Tvec = ...')
        logger.warning(numpy.std(Tvec, 0))

    return dict([( stats.items()[i][0], numpy.mean(stats.items()[i][1], 0) ) for i in xrange(0,len(stats))])


def calculate_Tz(mu0_try, t, delta):
    assert numpy.max(numpy.abs(t.diagonal())) == 0

    logger.debug(' mu0_try = %.9f', mu0_try)

    tt = t + mu0_try * numpy.identity(t.shape[0])
    soln = bcs_soln(tt, delta)
    stats = bcs_stats_average(soln['u'], soln['v'])

    logger.debug(' fdagf_try = %.9f', stats['fdagf'])

    return stats['Tvec'][2]


def calculate_pairing_matrix(u, v, norm):
    assert isinstance(u, numpy.ndarray)
    assert isinstance(v, numpy.ndarray)
    assert u.ndim == v.ndim == 2
    assert u.shape[0] == u.shape[1] == v.shape[0] == v.shape[1]

    # Warn if u is singular, but don't abort cuz sometimes this is OK..
    (signdet, logdet) = numpy.linalg.slogdet(u)
    if signdet * numpy.exp(logdet) == 0:
        logger.warning(' |det(u)| = %g', numpy.abs(numpy.exp(logdet)))

    phiMat = numpy.dot( numpy.linalg.inv(u.conjugate().transpose()), v.conjugate().transpose() )

#    logger.debug(' pairing matrix = ...')
#    logger.debug(phiMat)
    logger.debug(" pairing matrix's main diagonal = ...")
    logger.debug(phiMat.diagonal())

    # Check that phiMat is symmetric:
    if numpy.max(numpy.abs(phiMat - phiMat.transpose())) > 1e-8:
        logger.warning(' "symmetricness" of pairing matrix = %.9g', numpy.max(numpy.abs(phiMat - phiMat.transpose())))

    # Check against TI code (must first run *_Fourier.py with same parameters):
#    phiMat_TI = numpy.loadtxt('phiMat_TI.dat').view(complex)
#    logger.debug(' max|phiMat - phiMat_TI| = %.9g ', numpy.max(numpy.abs(phiMat - phiMat_TI)))

    phiNorm = norm_of_pairing_matrix(phiMat)
    logger.info(' norm of pairing matrix before normalization = %f', phiNorm)

    if norm is not None:
        phiMat = norm * (phiMat / phiNorm)
        phiNorm = norm_of_pairing_matrix(phiMat)
        logger.info(' norm of pairing matrix after normalization = %f', phiNorm)

    return phiMat


def norm_of_pairing_matrix(phiMat):
    return numpy.sqrt(numpy.sum(numpy.abs(phiMat)**2))


def plot_pairing_matrix(phiMat):
    assert isinstance(phiMat, numpy.ndarray)
    assert phiMat.ndim == 2
    assert phiMat.shape[0] == phiMat.shape[1]

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = numpy.arange(0, phiMat.shape[0], 1)
    Y = numpy.arange(0, phiMat.shape[0], 1)
    X, Y = numpy.meshgrid(X, Y)
    surf = ax.plot_surface(X, Y, numpy.abs(phiMat), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    #ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


def create_t_matrix_basic(lattice, boundary_conditions, t_value, neighbors_func):
    Nsites = len(lattice)

    t = numpy.zeros((Nsites, Nsites), dtype=complex)

    for site_r in lattice:
        ind_r = lattice.index(site_r)
        for site in neighbors_func(site_r):
            site_rp, phase = lattice.enforce_boundary(site, boundary_conditions)
            ind_rp = lattice.index(site_rp)
            t[ind_r, ind_rp] = t_value * phase

    assert numpy.max(numpy.abs(t - t.conjugate().transpose())) < 1e-12

    return t


def create_delta_matrix(lattice, boundary_conditions, delta_value, pairing_func, neighbors_func, channel):
    Nsites = len(lattice)

    delta = numpy.zeros((Nsites, Nsites), dtype=complex)

    for site_r in lattice:
        ind_r = lattice.index(site_r)
        for site in neighbors_func(site_r):
            site_rp, phase = lattice.enforce_boundary(site, boundary_conditions)
            ind_rp = lattice.index(site_rp)
            delta[ind_r, ind_rp] = delta_value * pairing_func(lattice, site_r, site) * phase

    if channel == 'singlet':
        assert numpy.max(numpy.abs(delta - delta.transpose())) < 1e-12
    elif channel == 'triplet':
        assert numpy.max(numpy.abs(delta + delta.transpose())) < 1e-12
    else:
        logger.warning(' unable to check delta for appropriate pairing channel (singlet or triplet).')

    return delta


##### Ansatz-specific stuff #####

def did_pairing_pattern(lattice, site1, site2):
    r1 = lattice.spatial_coordinates(site1)
    r2 = lattice.spatial_coordinates(site2)

    R = numpy.subtract(r2, r1)
    Rhat = R / numpy.linalg.norm(R)

    return (Rhat[0] + 1j*Rhat[1])**2


def dx2minusy2_pairing_pattern(lattice, site1, site2):
    r1 = lattice.spatial_coordinates(site1)
    r2 = lattice.spatial_coordinates(site2)

    R = numpy.subtract(r2, r1)
    Rhat = R / numpy.linalg.norm(R)

    return Rhat[0]**2 - Rhat[1]**2


def did_nn_bcs_theory(lattice, boundary_conditions, t1, delta1, delta0, mu0_start, mu0=None, norm=None):
    Nsites = len(lattice)

    t = create_t_matrix_basic(lattice, boundary_conditions, t1, lattice.nearest_neighbors)
    delta = create_delta_matrix(lattice, boundary_conditions, delta1, did_pairing_pattern, lattice.nearest_neighbors, 'singlet')

    assert numpy.max(numpy.abs(t.diagonal())) == 0
    assert numpy.max(numpy.abs(delta.diagonal())) == 0

    if mu0 is None:
        mu0_soln = optimize.fsolve(calculate_Tz, mu0_start, args=(t, delta), full_output=True, xtol=1e-06, epsfcn=0.1)
        logger.info(' mu0_soln = %s', mu0_soln)
        mu0 = mu0_soln[0][0]
#        mu0 = optimize.newton(calculate_Tz, mu0_start, tol=1e-06)

    tt = t + mu0 * numpy.identity(Nsites)
    delta = delta + delta0 * numpy.identity(Nsites)

    soln = bcs_soln(tt, delta)
    stats = bcs_stats_average(soln['u'], soln['v'])

    logger.info(' mu0 = %.9f', mu0)
    logger.info(' fdagf = %.9f', stats['fdagf'])
    logger.info(' Ttot = %.9f', stats['Ttot'])

    phi = calculate_pairing_matrix(soln['u'], soln['v'], norm)

    return {'pairing_matrix':phi, 'bcs_stats':stats, 'chemical_potential':mu0}


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    lattice = HexagonalLattice([6, 6])
    boundary_conditions = (periodic, antiperiodic)

    theta = 0.816814

    t1 = numpy.cos(theta)
    delta1 = numpy.sin(theta)

    bcs_theory = did_nn_bcs_theory(lattice, boundary_conditions, t1, delta1, 0, 0)

    plot_pairing_matrix(bcs_theory['pairing_matrix'])
