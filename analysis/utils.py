import numpy as np
from scipy.optimize import minimize, root
from time import perf_counter as time
from XCOM import mu_tot, mu_PE, mu_CS, mu_PP
import os
import re

def calcEnergyDeposited(filename, N):
    """Calculates the total energy deposited (and uncertainty) for filename"""
    E_deposited = 0
    sigma_deposited = 0
    with open(filename) as f:
        header = np.array(f.readline().split())
        idx = np.argmax(header == "E_deposited(MeV)")
        for line in f:
            energy = float(line.split()[idx])
            E_deposited += energy
            sigma_deposited += energy**2
    
    E_deposited = E_deposited / N
    sigma_deposited = np.sqrt(sigma_deposited) / N
    return E_deposited, sigma_deposited

def calcLookupTables(npy_file, *args):
    """Computes lookup tables from a directory of .dat files"""
    if len(args) == 0:
        lookup_alpha = np.load(npy_file, allow_pickle='TRUE').item()
        return lookup_alpha
    elif len(args) == 9:
        path, energies, thetaMax, R, E_dep, attenMat, zRange, compound_Z, compound_f = args
    else:
        path, energies, thetaMax, R, E_dep, attenMat, zRange = args
    t0 = time()
    print("Building lookup tables...")
    E_openBeam = {}
    sigma_openBeam = {}
    for filename in os.listdir(path):
        if filename.endswith(".dat"):
            E = re.search("E=(\d+\.?\d+|\d+)", filename)[1]
            lmbda = int(re.search("lmbda=(\d+\.?\d+|\d+)", filename)[1])
            if lmbda == 0 and E in energies:
                N = int(re.search("N=(\d+)", filename)[1])
                E_deposited, sigma_deposited = calcEnergyDeposited(path + filename, N)
                E_openBeam[E] = E_deposited
                sigma_openBeam[E] = sigma_deposited

    lookupTables = {E: {} for E in energies}
    for filename in os.listdir(path):
        if filename.endswith(".dat"):
            E = re.search("E=(\d+\.?\d+|\d+)", filename)[1]
            lmbda = int(re.search("lmbda=(\d+\.?\d+|\d+)", filename)[1])
            if lmbda > 0 and E in energies:
                N = int(re.search("N=(\d+)", filename)[1])
                Z = int(re.search("Z=(\d+)", filename)[1])
                if Z not in lookupTables[E].keys():
                    lookupTables[E][Z] = {}
                E_deposited, sigma_deposited = calcEnergyDeposited(path + filename, N)
                if E_deposited == 0:
                    print("zero energy deposited for (E, Z, lmbda) =", (E, Z, lmbda))
                    E_deposited = 1e-10
                    sigma_deposited = 1
                alpha = -np.log(E_deposited / E_openBeam[E])
                sigma = np.sqrt(sigma_deposited**2 / E_deposited**2 + sigma_openBeam[E]**2 / E_openBeam[E]**2)
                lookupTables[E][Z][lmbda] = (alpha, sigma)

    Z_arr = np.sort(np.array(list(lookupTables[energies[0]].keys())))
    lookup_alpha = {E: {} for E in energies}
    for Z in Z_arr:
        for E in energies:
            b = np.load("/Users/peter/Work/radiography/data/b%sMeV_10.npy" % E)
            lmbda_arr = np.sort(np.array(list(lookupTables[E][Z].keys())))
            alpha_arr, sigma_arr = np.array([lookupTables[E][Z][lmbda] for lmbda in lmbda_arr]).T

            if Z <= 100:
                lmbda_eff = np.array([calcLambdaEff(lmbda, thetaMax, Z, b, R, E_dep, attenMat, zRange) for lmbda in lmbda_arr])
            else:
                Z_arr = compound_Z[Z]
                f_arr = compound_f[Z]
                lmbda_eff = np.array([calcCompoundLambdaEff(lmbda, thetaMax, Z_arr, f_arr, b, R, E_dep, attenMat, zRange) for lmbda in lmbda_arr])

            lookup_alpha[E][Z] = (lmbda_eff, alpha_arr, sigma_arr)

    np.save(npy_file, lookup_alpha)
    print("Finished in %d seconds" % (time() - t0))
    return lookup_alpha

def extractFromTables(lookup_alpha, H, L):
    E0 = list(lookup_alpha.keys())[0]
    Z_keys = np.sort(np.array(list(lookup_alpha[E0].keys())))
    alpha_H_arr = []
    alpha_L_arr = []
    sigma_H_arr = []
    sigma_L_arr = []
    lmbda_H_arr = []
    lmbda_L_arr = []
    Z_arr = []
    Z_vals = Z_keys[Z_keys <= 100]
    Z_compound = Z_keys[Z_keys > 100]
    for Z in Z_vals:
        lmbda_H0, alpha_H0, sigma_H0 = lookup_alpha[H][Z]
        lmbda_L0, alpha_L0, sigma_L0 = lookup_alpha[L][Z]

        alpha_H_arr.append(alpha_H0)
        alpha_L_arr.append(alpha_L0)
        sigma_H_arr.append(sigma_H0)
        sigma_L_arr.append(sigma_L0)
        lmbda_H_arr.append(lmbda_H0)
        lmbda_L_arr.append(lmbda_L0)
        Z_arr.append(np.full(lmbda_H0.size, Z))

    alpha_H_arr = np.concatenate(alpha_H_arr)
    alpha_L_arr = np.concatenate(alpha_L_arr)
    sigma_H_arr = np.concatenate(sigma_H_arr)
    sigma_L_arr = np.concatenate(sigma_L_arr)
    lmbda_H_arr = np.concatenate(lmbda_H_arr)
    lmbda_L_arr = np.concatenate(lmbda_L_arr)
    Z_arr = np.concatenate(Z_arr)
    return Z_arr, Z_vals, Z_compound, alpha_H_arr, alpha_L_arr, sigma_H_arr, sigma_L_arr, lmbda_H_arr, lmbda_L_arr

def calcAttenMat_tot(E_g, zRange):
    """
    Returns a mass attenuation coefficient matrix, where element (i, j) is the mass
    attenuation coefficient of element with energy E_g[i] and atomic number zRange[j]
    """
    n = zRange.size
    attenMat = np.zeros((E_g.size, n))
    for i in range(n):
        Z = zRange[i]
        atten = mu_tot(E_g, Z)
        attenMat[:,i] = atten
    return attenMat

def calcAttenMat_PE(E_g, zRange):
    n = zRange.size
    attenMat_PE = np.zeros((E_g.size, n))
    for i in range(n):
        Z = zRange[i]
        atten = mu_PE(E_g, Z)
        attenMat_PE[:,i] = atten
    return attenMat_PE

def calcAttenMat_CS(E_g, zRange):
    n = zRange.size
    attenMat_CS = np.zeros((E_g.size, n))
    for i in range(n):
        Z = zRange[i]
        atten = mu_CS(E_g, Z)
        attenMat_CS[:,i] = atten
    return attenMat_CS

def calcAttenMat_PP(E_g, zRange):
    n = zRange.size
    attenMat_PP = np.zeros((E_g.size, n))
    for i in range(n):
        Z = zRange[i]
        atten = mu_PP(E_g, Z)
        attenMat_PP[:,i] = atten
    return attenMat_PP

def calcAlpha(lmbda, Z, b, R, E_dep, attenMat, zRange):
    """Calculate alpha for a given array of materials and thicknesses"""
    qb = R.T @ E_dep * b
    d = np.sum(qb)
    atten = attenMat[:,Z - zRange[0]]
    m0 = np.exp(-atten * lmbda)
    d0 = qb @ m0
    alpha = np.log(d / d0)
    return alpha

def calcCompoundAlpha(lmbda_arr, Z_arr, b, R, E_dep, attenMat, zRange):
    """Calculate alpha for a given compound material"""  
    qb = R.T @ E_dep * b
    d = np.sum(qb)
    atten = attenMat[:,Z_arr - zRange[0]]
    m0 = np.exp(-np.sum(atten * lmbda_arr, axis=1))
    d0 = qb @ m0
    alpha = np.log(d / d0)
    return alpha

def lookup(lmbda, Z, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange):
    """Calculate alpha and its derivatives for a given material and array of thicknesses"""
    q = R.T @ E_dep
    qb_H = q * b_H
    qb_L = q * b_L
    d_H = np.dot(q, b_H)
    d_L = np.dot(q, b_L)
    atten_H = attenMat_H[:,Z - zRange[0]]
    atten_L = attenMat_L[:,Z - zRange[0]]
    m0_H = np.exp(-np.outer(atten_H, lmbda))
    m1_H = -atten_H[:,None] * m0_H
    m2_H = -atten_H[:,None] * m1_H
    m0_L = np.exp(-np.outer(atten_L, lmbda))
    m1_L = -atten_L[:,None] * m0_L
    m2_L = -atten_L[:,None] * m1_L
    d_H0 = qb_H @ m0_H
    d_L0 = qb_L @ m0_L
    d_H1 = qb_H @ m1_H
    d_L1 = qb_L @ m1_L
    d_H2 = qb_H @ m2_H
    d_L2 = qb_L @ m2_L
    alpha_H0 = np.log(d_H / d_H0)
    alpha_L0 = np.log(d_L / d_L0)
    alpha_H1 = -d_H1 / d_H0
    alpha_L1 = -d_L1 / d_L0
    alpha_H2 = (d_H1**2 - d_H0 * d_H2) / d_H0**2
    alpha_L2 = (d_L1**2 - d_L0 * d_L2) / d_L0**2
    return alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2

def fitToTheory(alpha_H, alpha_L, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange):
    """For a list of alpha_H, alpha_L values, finds the theoretical Z which best reproduces the list"""
    b = zRange.size
    loss_arr = np.zeros(b)
    lmbda = np.ones(alpha_H.size)
    for i in range(b):
        Z = zRange[i]
        if i == 0:
            nsteps = 6
        elif i == 1:
            nsteps = 4
        else:
            nsteps = 2  
        for _ in range(nsteps):
            alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookup(lmbda, Z, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange)
            diff_H = alpha_H0 - alpha_H
            diff_L = alpha_L0 - alpha_L

            grad = diff_H * alpha_H1 + diff_L * alpha_L1
            hess = diff_H * alpha_H2 + alpha_H1**2 + diff_L * alpha_L2 + alpha_L1**2
            lmbda = lmbda - grad / hess

        alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookup(lmbda, Z, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange)
        loss_vec = (alpha_H0 - alpha_H)**2 + (alpha_L0 - alpha_L)**2
        loss_arr[i] = np.mean(loss_vec)
    
    loss_pad = np.pad(loss_arr, 1, constant_values = np.inf)
    left = loss_pad[:-2]
    right = loss_pad[2:]
    optima = (left > loss_arr) & (right > loss_arr)
    thresh = loss_arr < 3*np.min(loss_arr)
    idx = np.argwhere(optima & thresh).flatten()
    while idx.size > 2:
        idx = np.delete(idx, np.argmax(loss_arr[idx]))
    return zRange[idx], loss_arr[idx]

### Test debug

# def lookupFrac(lmbda, Z, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange):
#     """Calculate alpha and its derivatives for a given material and array of thicknesses"""
#     q = R.T @ E_dep
#     qb_H = q * b_H
#     qb_L = q * b_L
#     d_H = np.dot(q, b_H)
#     d_L = np.dot(q, b_L)
#     Z_floor = np.floor(Z).astype(int)
#     Z_ceil = np.ceil(Z).astype(int)
#     f = Z_ceil - Z
#     atten_H_floor = attenMat_H[:,Z_floor - zRange[0]]
#     atten_H_ceil = attenMat_H[:,Z_ceil - zRange[0]]
#     atten_H = f * atten_H_floor + (1 - f) * atten_H_ceil
#     atten_L_floor = attenMat_L[:,Z_floor - zRange[0]]
#     atten_L_ceil = attenMat_L[:,Z_ceil - zRange[0]]
#     atten_L = f * atten_L_floor + (1 - f) * atten_L_ceil
#     m0_H = np.exp(-np.outer(atten_H, lmbda))
#     m1_H = -atten_H[:,None] * m0_H
#     m2_H = -atten_H[:,None] * m1_H
#     m0_L = np.exp(-np.outer(atten_L, lmbda))
#     m1_L = -atten_L[:,None] * m0_L
#     m2_L = -atten_L[:,None] * m1_L
#     d_H0 = qb_H @ m0_H
#     d_L0 = qb_L @ m0_L
#     d_H1 = qb_H @ m1_H
#     d_L1 = qb_L @ m1_L
#     d_H2 = qb_H @ m2_H
#     d_L2 = qb_L @ m2_L
#     alpha_H0 = np.log(d_H / d_H0)
#     alpha_L0 = np.log(d_L / d_L0)
#     alpha_H1 = -d_H1 / d_H0
#     alpha_L1 = -d_L1 / d_L0
#     alpha_H2 = (d_H1**2 - d_H0 * d_H2) / d_H0**2
#     alpha_L2 = (d_L1**2 - d_L0 * d_L2) / d_L0**2
#     return alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2

# def fitToTheory(alpha_H, alpha_L, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange):
#     """For a list of alpha_H, alpha_L values, finds the theoretical Z which best reproduces the list"""
#     zRangeFrac = np.linspace(zRange[0], zRange[-1], 10*(zRange.size-1) + 1)
#     b = zRangeFrac.size
#     loss_arr = np.zeros(b)
#     lmbda = np.ones(alpha_H.size)
#     for i in range(b):
#         Z = zRangeFrac[i]
#         if i == 0:
#             nsteps = 4
#         else:
#             nsteps = 1
#         for _ in range(nsteps):
#             alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookupFrac(lmbda, Z, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange)
#             diff_H = alpha_H0 - alpha_H
#             diff_L = alpha_L0 - alpha_L

#             grad = diff_H * alpha_H1 + diff_L * alpha_L1
#             hess = diff_H * alpha_H2 + alpha_H1**2 + diff_L * alpha_L2 + alpha_L1**2
#             lmbda = lmbda - grad / hess

#         alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookupFrac(lmbda, Z, b_H, b_L, R, E_dep, attenMat_H, attenMat_L, zRange)
#         loss_vec = (alpha_H0 - alpha_H)**2 + (alpha_L0 - alpha_L)**2
#         loss_arr[i] = np.mean(loss_vec)

#     idx = np.argmin(loss_arr)
#     return zRangeFrac[idx], loss_arr[idx]

# def calcAlphaFrac(lmbda, Z, b, R, E_dep, attenMat, zRange):
#     """Calculate alpha for a given array of materials and thicknesses"""
#     qb = R.T @ E_dep * b
#     d = np.sum(qb)
#     Z_floor = np.floor(Z).astype(int)
#     Z_ceil = np.ceil(Z).astype(int)
#     f = Z_ceil - Z
#     atten_floor = attenMat[:,Z_floor - zRange[0]]
#     atten_ceil = attenMat[:,Z_ceil - zRange[0]]
#     atten = f * atten_floor + (1 - f) * atten_ceil
#     m0 = np.exp(-atten * lmbda)
#     d0 = qb @ m0
#     alpha = np.log(d / d0)
#     return alpha

###

def calcLambdaEff(lmbda0, thetaMax, Z, b, R, E_dep, attenMat, zRange):
    """Finds the effective lambda which approximates the entire target"""
    def calcLogD0(lmbda):
        m0 = np.exp(-atten * lmbda)
        d0 = np.dot(q, m0 * b)
        logD0 = np.log(d0)
        return logD0
        
    atten = attenMat[:,Z - zRange[0]]
    thetaRange = np.linspace(0, thetaMax, 1001)
    dtheta = thetaRange[1] - thetaRange[0]
    theta = 0.5*(thetaRange[1:] + thetaRange[:-1])
    lmbda = lmbda0 / np.cos(theta)
    m = np.exp(-np.outer(atten, lmbda))
    m0 = dtheta * np.sum(m, axis=1) / thetaMax
    q = R.T @ E_dep
    d0 = np.dot(q, m0 * b)
    logD0 = np.log(d0)
    
    sol = root(lambda lmbdaEff: calcLogD0(lmbdaEff) - logD0, x0=lmbda0)
    assert sol.success
    return sol.x[0]

def calcCompoundLambdaEff(lmbda0, thetaMax, Z_arr, f_arr, b, R, E_dep, attenMat, zRange):
    """Finds the effective lambda which approximates the entire compound target"""
    def calcLogD0(lmbda):
        m0 = np.exp(-np.sum(atten * (lmbda * f_arr), axis=1))
        d0 = np.dot(q, m0 * b)
        logD0 = np.log(d0)
        return logD0
    
    atten = attenMat[:,Z_arr - zRange[0]]
    thetaRange = np.linspace(0, thetaMax, 1001)
    dtheta = thetaRange[1] - thetaRange[0]
    theta = 0.5*(thetaRange[1:] + thetaRange[:-1])
    A = np.sum(atten * lmbda0 * f_arr, axis=1)
    m = np.exp(-np.outer(A, 1/np.cos(theta)))
    m0 = dtheta * np.sum(m, axis=1) / thetaMax
    q = R.T @ E_dep
    d0 = np.dot(q, m0 * b)
    logD0 = np.log(d0)
    
    sol = root(lambda lmbdaEff: calcLogD0(lmbdaEff) - logD0, x0=lmbda0)
    assert sol.success
    return sol.x[0]

def calcBias(alpha, sigma, lmbda, Z, b, R, E_dep, attenMat, zRange):
    """Finds the bias term such that chi-squared equals one"""
    def calcChiSquared(bias):
        alpha0 = calcAlpha(lmbda, Z, b, R, E_dep, attenMat, zRange)
        chi2_vec = (alpha0 - alpha)**2 / (sigma**2 + bias**2)
        return np.mean(chi2_vec)

    sol = root(lambda bias: calcChiSquared(bias)-1, x0=0.1)
    assert sol.success
    return np.abs(sol.x)

def fitSemiempirical(alpha, lmbda, Z, b, R, E_dep, attenMat_tot, attenMat_PE, attenMat_CS, attenMat_PP, zRange):
    """For a list of alpha values, finds best 'a', 'b', and 'c' coefficient to reproduce 
    the given lambda and Z values"""
    def calcLoss(x):
        attenMat = attenMat_tot + (x[0]-1)*attenMat_PE + (x[1]-1)*attenMat_CS + (x[2]-1)*attenMat_PP
        alpha0 = calcAlpha(lmbda, Z, b, R, E_dep, attenMat, zRange)
        loss_vec = (alpha0 - alpha)**2
        return np.mean(loss_vec)
    
    res = minimize(calcLoss, x0=(1, 1, 1))
    assert res.success
    a0, b0, c0 = res.x
    loss = res.fun
    print("Minimum found at a = %.4f, b=%.4f, c=%.4f with a loss of %.3e" % (a0, b0, c0, loss))
    return a0, b0, c0
