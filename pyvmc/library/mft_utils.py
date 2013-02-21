#!/usr/bin/env python

from pyvmc.core import Lattice, LatticeSite, HexagonalLattice, periodic, antiperiodic
from pyvmc.utils import is_square_matrix

import numpy
from numpy import pi
numpy.set_printoptions(precision=15, suppress=True, threshold=10000, linewidth=300)

from scipy import optimize

import logging
logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt


##### General stuff #####

def free_fermion_soln(t):
    assert is_square_matrix(t)

    Nsites = t.shape[0]

    # t_{ij} = ttilde_{ij} = t_{ij} + mu_i * krondelta_{ij} from notes (FIXME: take t and mu as separate parameters?)
    H = -t

    assert numpy.max(numpy.abs(H - H.conjugate().transpose())) < 1e-12

    eigsys = numpy.linalg.eigh(H)

    idx = eigsys[0].argsort()
    evals = eigsys[0][idx]
    evecs = eigsys[1][:,idx]

    logger.debug('Eigenvalues = ...')
    logger.debug(evals)

    return {'evecs':evecs, 'evals':evals}


def calculate_free_fermion_pairing_matrix(evecs, evals, Nparticles, norm):
    assert is_square_matrix(evecs)

    Nsites = evecs.shape[0]
    assert Nsites == len(evals)

    M = Nparticles/2  # ( = N_up = N_down )
    assert M <= Nsites

    phiMat = numpy.zeros((Nsites, Nsites), dtype=complex)
    for ind_r in range(0, Nsites):
        for ind_rp in range(0, Nsites):
            phiMat[ind_r, ind_rp] = numpy.dot( evecs[ind_r, 0:M], evecs[ind_rp, 0:M] )

#    logger.debug('Pairing matrix = ...')
#    logger.debug(phiMat)
#    logger.debug("Pairing matrix's main diagonal = ...")
#    logger.debug(phiMat.diagonal())

    # Check that phiMat is symmetric:
    if numpy.max(numpy.abs(phiMat - phiMat.transpose())) > 1e-8:
        raise RuntimeError('"Symmetricness" of pairing matrix = %.9g', numpy.max(numpy.abs(phiMat - phiMat.transpose())))

    phiNorm = norm_of_pairing_matrix(phiMat)
    logger.info('Norm of pairing matrix before normalization = %f', phiNorm)

    if norm is not None:
        phiMat = norm * (phiMat / phiNorm)
        phiNorm = norm_of_pairing_matrix(phiMat)
        logger.info('Norm of pairing matrix after normalization = %f', phiNorm)

    return phiMat


def bcs_soln(t, delta):
    # both t and delta should include on-site terms, e.g., t_{ij} = ttilde_{ij} = t_{ij} + mu_i * krondelta_{ij} from notes
    assert is_square_matrix(t)
    assert is_square_matrix(delta)
    assert t.shape[0] == delta.shape[0]

    Nsites = t.shape[0]

    H = numpy.zeros((Nsites, Nsites, 2, 2), dtype=complex)
    for ind_r in range(0, Nsites):
        for ind_rp in range(0, Nsites):
            H[ind_r, ind_rp] = -numpy.array([[t[ind_r, ind_rp], delta[ind_r, ind_rp]],
                                             [delta[ind_r, ind_rp].conjugate(), -t[ind_r, ind_rp].conjugate()]])

    Hflat = numpy.zeros((2*Nsites, 2*Nsites), dtype=complex)
    for i in range(0, 2*Nsites):
        for j in range(0, 2*Nsites):
            Hflat[i,j] = H[i//2, j//2, numpy.mod(i,2), numpy.mod(j,2)]

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

    logger.debug('Positive eigenvalues = ...')
    logger.debug(evalsPos)
#    logger.debug('Negative eigenvalues = ...')
#    logger.debug(evalsNeg[::-1])

    # Check that positive eigenvalues are positive and that eigenvalues have come in plus-minus pairs:
    # FIXME:  probably shouldn't be so strict, since this gets called when trying to tune the ansatz,
    # and we could potentially run into degenerate mean field solutions (i.e., those with zero modes) in that process;
    # if we don't assert evalsPos > 0, then we treat zero modes as positive and don't populate them in the the BCS solution,
    # which I think is OK; for now, we'll just warn instead of assert
#    assert numpy.min(evalsPos) > 1e-12
    if not numpy.min(evalsPos) > 1e-12:
        logger.warning('Some of evalsPos are not greater than 1e-12')
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
    assert is_square_matrix(u)
    assert is_square_matrix(v)
    assert u.shape[0] == v.shape[0]

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
        logger.warning('This BCS state is not translationally invariant!  Standard deviation of Tvec = ...')
        logger.warning(numpy.std(stats['Tvec'], 0))

    return {k: numpy.mean(v, 0) for k, v in stats.items()}


def calculate_Tz(mu0_try, t_offsite, delta):
    logger.debug('mu0_try = %.9f', mu0_try)

    t = add_onsite_terms(t_offsite, mu0_try)
    soln = bcs_soln(t, delta)
    stats = bcs_stats_average(soln['u'], soln['v'])

    logger.debug('fdagf_try = %.9f', stats['fdagf'])

    return stats['Tvec'][2]


def calculate_pairing_matrix(u, v, norm):
    assert is_square_matrix(u)
    assert is_square_matrix(v)
    assert u.shape[0] == v.shape[0]

    # Warn if u is singular, but don't abort cuz sometimes this is OK..
    (signdet, logdet) = numpy.linalg.slogdet(u)
    logger.info('|det(u)| = %g', numpy.abs(numpy.exp(logdet)))  # actually, just log it for now..

    phiMat = numpy.dot( numpy.linalg.inv(u.conjugate().transpose()), v.conjugate().transpose() )

#    logger.debug('Pairing matrix = ...')
#    logger.debug(phiMat)
#    logger.debug("Pairing matrix's main diagonal = ...")
#    logger.debug(phiMat.diagonal())

    # Check that phiMat is symmetric:
    if numpy.max(numpy.abs(phiMat - phiMat.transpose())) > 1e-8:
        raise RuntimeError('"Symmetricness" of pairing matrix = %.9g', numpy.max(numpy.abs(phiMat - phiMat.transpose())))

    # Check against TI code (must first run *_Fourier.py with same parameters):
#    phiMat_TI = numpy.loadtxt('phiMat_TI.dat').view(complex)
#    logger.debug('max|phiMat - phiMat_TI| = %.9g ', numpy.max(numpy.abs(phiMat - phiMat_TI)))

    phiNorm = norm_of_pairing_matrix(phiMat)
    logger.info('Norm of pairing matrix before normalization = %f', phiNorm)

    if norm is not None:
        phiMat = norm * (phiMat / phiNorm)
        phiNorm = norm_of_pairing_matrix(phiMat)
        logger.info('Norm of pairing matrix after normalization = %f', phiNorm)

    return phiMat


def norm_of_pairing_matrix(phiMat):
    return numpy.sqrt(numpy.sum(numpy.abs(phiMat)**2))


def plot_pairing_matrix(phiMat, callshow=True):
    assert is_square_matrix(phiMat)

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

    if callshow:
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

    check_t_matrix(t)  # FIXME: not sure if we always want to check this here, or outside the function

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

    check_delta_matrix(delta, channel)  # FIXME: not sure if we always want to check this here, or outside the function

    return delta


def add_onsite_terms(offsite, onsite):
    assert is_square_matrix(offsite)

    check_offsite_for_onsite(offsite)  # since we're singling out the on-site terms with this function

    Nsites = offsite.shape[0]

    # total_{ij} = offsite_{ij} + onsite_i * krondelta_{ij}
    if isinstance(onsite, numpy.ndarray):
        assert onsite.ndim == 1
        assert onsite.shape[0] == 1 or onsite.shape[0] == Nsites
        if onsite.shape[0] == 1:
            total = offsite + onsite * numpy.identity(Nsites)
        elif onsite.shape[0] == Nsites:
            total = offsite + numpy.diag(onsite)
    else:
        total = offsite + onsite * numpy.identity(Nsites)

    return total


def check_t_matrix(t):
    assert numpy.max(numpy.abs(t - t.conjugate().transpose())) < 1e-12


def check_delta_matrix(delta, channel):
    if channel == 'singlet':
        assert numpy.max(numpy.abs(delta - delta.transpose())) < 1e-12
    elif channel == 'triplet':
        assert numpy.max(numpy.abs(delta + delta.transpose())) < 1e-12
    else:
        logger.warning('Unable to check delta for appropriate pairing channel (singlet or triplet).')


def check_offsite_for_onsite(offsite):
    assert numpy.max(numpy.abs(offsite.diagonal())) == 0



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


# hf is for 'half-filling' since the function currently solves for the root of Tz: Tz=0 (fdagf=1)
# FIXME: have desired average fdagf or Tvec be a parameter
def did_hf_bcs_theory(lattice, boundary_conditions, t1, delta1, mu0_start=0, mu0=None, delta0=0, t2=0, delta2=0, t3=0, delta3=0):

    delta_offsite = \
        create_delta_matrix(lattice, boundary_conditions, delta1, did_pairing_pattern, lattice.nearest_neighbors, 'singlet') + \
        create_delta_matrix(lattice, boundary_conditions, delta2, did_pairing_pattern, lattice.second_nearest_neighbors, 'singlet') + \
        create_delta_matrix(lattice, boundary_conditions, delta3, did_pairing_pattern, lattice.third_nearest_neighbors, 'singlet')
    delta = add_onsite_terms(delta_offsite, delta0)
    check_delta_matrix(delta, 'singlet')  # this will never fail for singlet pairing (create_delta_matrix also checks delta)

    t_offsite = \
        create_t_matrix_basic(lattice, boundary_conditions, t1, lattice.nearest_neighbors) + \
        create_t_matrix_basic(lattice, boundary_conditions, t2, lattice.second_nearest_neighbors) + \
        create_t_matrix_basic(lattice, boundary_conditions, t3, lattice.third_nearest_neighbors)
    check_offsite_for_onsite(t_offsite)

    if mu0 is None:
        mu0_soln = optimize.fsolve(calculate_Tz, mu0_start, args=(t_offsite, delta_offsite), full_output=True, xtol=1e-06, epsfcn=0.1)
        logger.info('mu0_soln = %s', mu0_soln)
        mu0 = mu0_soln[0][0]
#        mu0 = optimize.newton(calculate_Tz, mu0_start, tol=1e-06)

    t = add_onsite_terms(t_offsite, mu0)
    check_t_matrix(t)  # this would only fail if mu0 were complex or something

    soln = bcs_soln(t, delta)
    stats = bcs_stats_average(soln['u'], soln['v'])

    logger.info('mu0 = %.9f', mu0)
    logger.info('fdagf = %.9f', stats['fdagf'])
    logger.info('Ttot = %.9f', stats['Ttot'])

    phi = calculate_pairing_matrix(soln['u'], soln['v'], norm=len(lattice))

    return {'pairing_matrix':phi, 'chemical_potential':mu0, 'bcs_stats':stats}


# hf is for 'half-filling' since the function currently solves for the root of Tz: Tz=0 (fdagf=1)
# FIXME: have desired average fdagf or Tvec be a parameter
def dx2minusy2_hf_bcs_theory(lattice, boundary_conditions, t1, delta1, mu0_start=0, mu0=None, delta0=0, t2=0, delta2=0, t3=0, delta3=0):

    delta_offsite = \
        create_delta_matrix(lattice, boundary_conditions, delta1, dx2minusy2_pairing_pattern, lattice.nearest_neighbors, 'singlet') + \
        create_delta_matrix(lattice, boundary_conditions, delta2, dx2minusy2_pairing_pattern, lattice.second_nearest_neighbors, 'singlet') + \
        create_delta_matrix(lattice, boundary_conditions, delta3, dx2minusy2_pairing_pattern, lattice.third_nearest_neighbors, 'singlet')
    delta = add_onsite_terms(delta_offsite, delta0)
    check_delta_matrix(delta, 'singlet')  # this will never fail for singlet pairing (create_delta_matrix also checks delta)

    t_offsite = \
        create_t_matrix_basic(lattice, boundary_conditions, t1, lattice.nearest_neighbors) + \
        create_t_matrix_basic(lattice, boundary_conditions, t2, lattice.second_nearest_neighbors) + \
        create_t_matrix_basic(lattice, boundary_conditions, t3, lattice.third_nearest_neighbors)
    check_offsite_for_onsite(t_offsite)

    if mu0 is None:
        mu0_soln = optimize.fsolve(calculate_Tz, mu0_start, args=(t_offsite, delta_offsite), full_output=True, xtol=1e-06, epsfcn=0.1)
        logger.info('mu0_soln = %s', mu0_soln)
        mu0 = mu0_soln[0][0]
#        mu0 = optimize.newton(calculate_Tz, mu0_start, tol=1e-06)

    t = add_onsite_terms(t_offsite, mu0)
    check_t_matrix(t)  # this would only fail if mu0 were complex or something

    soln = bcs_soln(t, delta)
    stats = bcs_stats_average(soln['u'], soln['v'])

    logger.info('mu0 = %.9f', mu0)
    logger.info('fdagf = %.9f', stats['fdagf'])
    logger.info('Ttot = %.9f', stats['Ttot'])

    phi = calculate_pairing_matrix(soln['u'], soln['v'], norm=len(lattice))

    return {'pairing_matrix':phi, 'chemical_potential':mu0, 'bcs_stats':stats}


# FIXME: call signature should look like next line, but I guess could have it take a bunch of garbage arguments to fit into framework of traditional BCS states..
#def sbm_bcs_theory(lattice, boundary_conditions, t1=1, t2=0, t3=0):
def sbm_bcs_theory(lattice, boundary_conditions, t1=1, delta1=0, mu0_start=0, mu0=None, delta0=0, t2=0, delta2=0, t3=0, delta3=0):
    assert len(lattice) % 2 == 0

    t_offsite = \
        create_t_matrix_basic(lattice, boundary_conditions, t1, lattice.nearest_neighbors) + \
        create_t_matrix_basic(lattice, boundary_conditions, t2, lattice.second_nearest_neighbors) + \
        create_t_matrix_basic(lattice, boundary_conditions, t3, lattice.third_nearest_neighbors)
    check_offsite_for_onsite(t_offsite)

    t = t_offsite
    check_t_matrix(t)

    soln = free_fermion_soln(t)
    stats = {
        'fdagf': 1.0,
        'fdownfup': 0,
        'Tvec': numpy.array([0, 0, 0]),
        'Ttot': 0,
    }

    M = len(lattice)/2  # each species is at half-filling --> spin wave function
    mu0 = soln['evals'][M-1]  # energies are returned sorted in increasing order
#    print M, soln['evals'][(M-1)-1], soln['evals'][(M-1)+1]

    logger.info('mu0 = %.9f', mu0)
    logger.info('fdagf = %.9f', stats['fdagf'])
    logger.info('Ttot = %.9f', stats['Ttot'])

    phi = calculate_free_fermion_pairing_matrix(soln['evecs'], soln['evals'], Nparticles=len(lattice), norm=len(lattice))

    return {'pairing_matrix':phi, 'chemical_potential':mu0, 'bcs_stats':stats}



if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        )

    lattice = HexagonalLattice([12, 12])
    boundary_conditions = (periodic, antiperiodic)

#    alpha = pi/4

#    bcs_theory = dx2minusy2_hf_bcs_theory(lattice, boundary_conditions, t1=numpy.cos(alpha), delta1=numpy.sin(alpha), mu0_start=0)
#    bcs_theory = dx2minusy2_hf_bcs_theory(lattice, boundary_conditions, t1=1, delta1=5, mu0_start=0, mu0=None)

#    bcs_theory = did_hf_bcs_theory(lattice, boundary_conditions, t1=numpy.cos(alpha), delta1=numpy.sin(alpha), mu0=0.1)
#    bcs_theory = did_hf_bcs_theory(lattice, boundary_conditions, t1=1, delta1=5)
#    bcs_theory = did_hf_bcs_theory(lattice, boundary_conditions, t1=0, delta1=0.99, mu0=0, delta2=-1.01)

#    bcs_theory = did_hf_bcs_theory(lattice, boundary_conditions, t1=1, delta1=0, mu0_start=0.8, mu0=None)
    bcs_theory = sbm_bcs_theory(lattice, boundary_conditions, t1=1)

    plot_pairing_matrix(bcs_theory['pairing_matrix'])
