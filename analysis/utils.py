import numpy as np
from scipy.optimize import minimize, root
from XCOM import mu_tot, mu_PE, mu_CS, mu_PP, mu_en
import os

def calcLookupTables(data_dir):
    """Computes lookup tables from a directory of .npy files"""
    path = "/Users/peter/Work/semiempirical_transparency/data/"
    E_openBeam = {}
    var_openBeam = {}
    for filename in os.listdir(path + "open_beam"):
        if filename.endswith(".npy"):
            data = np.load(path + "open_beam/" + filename, allow_pickle=True).item()
            if data["lambda"] == 0:
                E = data["E"]
                E_openBeam[E] = data["E_deposited"]
                var_openBeam[E] = data["var_deposited"]

    energies = list(E_openBeam.keys())
    lookupTables = {E: {} for E in energies}
    for filename in os.listdir(path + data_dir):
        if filename.endswith(".npy"):
            data = np.load(path + data_dir + "/" + filename, allow_pickle=True).item()
            E = data["E"]
            E_deposited = data["E_deposited"]
            var_deposited = data["var_deposited"]
            lmbda_eff = data["lambda_eff"]
            Z = data["Z"]
            alpha = -np.log(E_deposited / E_openBeam[E])
            sigma = np.sqrt(var_deposited / E_deposited**2 + var_openBeam[E] / E_openBeam[E]**2)
            if Z not in lookupTables[E].keys():
                lookupTables[E][Z] = {}
            lookupTables[E][Z][lmbda_eff] = (alpha, sigma)

    Z_arr = np.sort(np.array(list(lookupTables[energies[0]].keys())))
    lookup_alpha = {E: {} for E in energies}
    for E in energies:
        for Z in Z_arr:
            lmbda_arr = np.sort(np.array(list(lookupTables[E][Z].keys())))
            alpha_arr, sigma_arr = np.array([lookupTables[E][Z][lmbda] for lmbda in lmbda_arr]).T
            lookup_alpha[E][Z] = (lmbda_arr, alpha_arr, sigma_arr)
    
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

def calcMu_tot(E, Z_range):
    """
    Returns a mass muuation coefficient matrix, where element (i, j) is the mass
    attenuation coefficient of element with energy E[i] and atomic number Z_range[j]
    """
    mu_mat = np.zeros((E.size, Z_range.size))
    for i in range(Z_range.size):
        Z = Z_range[i]
        mu_mat[:,i] = mu_tot(E, Z)
    return mu_mat

def calcMu_PE(E, Z_range):
    mu_mat_PE = np.zeros((E.size, Z_range.size))
    for i in range(Z_range.size):
        mu_mat_PE[:,i] = mu_PE(E, Z_range[i])
    return mu_mat_PE

def calcMu_CS(E, Z_range):
    mu_mat_CS = np.zeros((E.size, Z_range.size))
    for i in range(Z_range.size):
        mu_mat_CS[:,i] = mu_CS(E, Z_range[i])
    return mu_mat_CS

def calcMu_PP(E, Z_range):
    mu_mat_PP = np.zeros((E.size, Z_range.size))
    for i in range(Z_range.size):
        mu_mat_PP[:,i] = mu_PP(E, Z_range[i])
    return mu_mat_PP

def calcAlpha(lmbda, Z, phi, D, mu_mat, Z_range):
    """Calculate alpha for a given array of materials and thicknesses"""
    D_phi = D * phi
    d = np.sum(D_phi)
    mu = mu_mat[:,Z - Z_range[0]]
    m0 = np.exp(-mu * lmbda)
    d0 = D_phi @ m0
    alpha = np.log(d / d0)
    return alpha

def calcCompoundAlpha(lmbda_arr, Z_arr, phi, D, mu_mat, Z_range):
    """Calculate alpha for a given compound material"""  
    D_phi = D * phi
    d = np.sum(D_phi)
    mu = mu_mat[:,Z_arr - Z_range[0]]
    m0 = np.exp(-np.sum(mu * lmbda_arr, axis=1))
    d0 = D_phi @ m0
    alpha = np.log(d / d0)
    return alpha

def lookup(lmbda, Z, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range):
    """Calculate alpha and its derivatives for a given material and array of thicknesses"""
    D_phi_H = D * phi_H
    D_phi_L = D * phi_L
    d_H = np.dot(D, phi_H)
    d_L = np.dot(D, phi_L)
    mu_H = mu_mat_H[:,Z - Z_range[0]]
    mu_L = mu_mat_L[:,Z - Z_range[0]]
    m0_H = np.exp(-np.outer(mu_H, lmbda))
    m1_H = -mu_H[:,None] * m0_H
    m2_H = -mu_H[:,None] * m1_H
    m0_L = np.exp(-np.outer(mu_L, lmbda))
    m1_L = -mu_L[:,None] * m0_L
    m2_L = -mu_L[:,None] * m1_L
    d_H0 = D_phi_H @ m0_H
    d_L0 = D_phi_L @ m0_L
    d_H1 = D_phi_H @ m1_H
    d_L1 = D_phi_L @ m1_L
    d_H2 = D_phi_H @ m2_H
    d_L2 = D_phi_L @ m2_L
    alpha_H0 = np.log(d_H / d_H0)
    alpha_L0 = np.log(d_L / d_L0)
    alpha_H1 = -d_H1 / d_H0
    alpha_L1 = -d_L1 / d_L0
    alpha_H2 = (d_H1**2 - d_H0 * d_H2) / d_H0**2
    alpha_L2 = (d_L1**2 - d_L0 * d_L2) / d_L0**2
    return alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2

def fitToTheory(alpha_H, alpha_L, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range):
    """For a list of alpha_H, alpha_L values, finds the theoretical Z which best reproduces the list"""
    b = Z_range.size
    loss_arr = np.zeros(b)
    lmbda = np.ones(alpha_H.size)
    for i in range(b):
        Z = Z_range[i]
        if i == 0:
            nsteps = 6
        elif i == 1:
            nsteps = 4
        else:
            nsteps = 2  
        for _ in range(nsteps):
            alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookup(lmbda, Z, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range)
            diff_H = alpha_H0 - alpha_H
            diff_L = alpha_L0 - alpha_L

            grad = diff_H * alpha_H1 + diff_L * alpha_L1
            hess = diff_H * alpha_H2 + alpha_H1**2 + diff_L * alpha_L2 + alpha_L1**2
            lmbda = lmbda - grad / hess

        alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookup(lmbda, Z, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range)
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
    return Z_range[idx], loss_arr[idx]

### non-integer atomic numbers

# def lookupFrac(lmbda, Z, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range):
#     """Calculate alpha and its derivatives for a given material and array of thicknesses"""
#     D_phi_H = D * phi_H
#     D_phi_L = D * phi_L
#     d_H = np.dot(D, phi_H)
#     d_L = np.dot(D, phi_L)
#     Z_floor = np.floor(Z).astype(int)
#     Z_ceil = np.ceil(Z).astype(int)
#     f = Z_ceil - Z
#     mu_H_floor = mu_mat_H[:,Z_floor - Z_range[0]]
#     mu_H_ceil = mu_mat_H[:,Z_ceil - Z_range[0]]
#     mu_H = f * mu_H_floor + (1 - f) * mu_H_ceil
#     mu_L_floor = mu_mat_L[:,Z_floor - Z_range[0]]
#     mu_L_ceil = mu_mat_L[:,Z_ceil - Z_range[0]]
#     mu_L = f * mu_L_floor + (1 - f) * mu_L_ceil
#     m0_H = np.exp(-np.outer(mu_H, lmbda))
#     m1_H = -mu_H[:,None] * m0_H
#     m2_H = -mu_H[:,None] * m1_H
#     m0_L = np.exp(-np.outer(mu_L, lmbda))
#     m1_L = -mu_L[:,None] * m0_L
#     m2_L = -mu_L[:,None] * m1_L
#     d_H0 = D_phi_H @ m0_H
#     d_L0 = D_phi_L @ m0_L
#     d_H1 = D_phi_H @ m1_H
#     d_L1 = D_phi_L @ m1_L
#     d_H2 = D_phi_H @ m2_H
#     d_L2 = D_phi_L @ m2_L
#     alpha_H0 = np.log(d_H / d_H0)
#     alpha_L0 = np.log(d_L / d_L0)
#     alpha_H1 = -d_H1 / d_H0
#     alpha_L1 = -d_L1 / d_L0
#     alpha_H2 = (d_H1**2 - d_H0 * d_H2) / d_H0**2
#     alpha_L2 = (d_L1**2 - d_L0 * d_L2) / d_L0**2
#     return alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2

# def calcAlphaFrac(lmbda, Z, phi, D, mu_mat, Z_range):
#     """Calculate alpha for a given array of materials and thicknesses"""
#     D_phi = D * phi
#     d = np.sum(D_phi)
#     Z_floor = np.floor(Z).astype(int)
#     Z_ceil = np.ceil(Z).astype(int)
#     f = Z_ceil - Z
#     mu_floor = mu_mat[:,Z_floor - Z_range[0]]
#     mu_ceil = mu_mat[:,Z_ceil - Z_range[0]]
#     mu = f * mu_floor + (1 - f) * mu_ceil
#     m0 = np.exp(-mu * lmbda)
#     d0 = D_phi @ m0
#     alpha = np.log(d / d0)
#     return alpha

# def fitToTheoryFrac(alpha_H, alpha_L, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range):
#     """For a list of alpha_H, alpha_L values, finds the theoretical Z which best reproduces the list"""
#     Z_rangeFrac = np.linspace(Z_range[0], Z_range[-1], 10*(Z_range.size-1) + 1)
#     b = Z_rangeFrac.size
#     loss_arr = np.zeros(b)
#     lmbda = np.ones(alpha_H.size)
#     for i in range(b):
#         Z = Z_rangeFrac[i]
#         if i == 0:
#             nsteps = 4
#         else:
#             nsteps = 1
#         for _ in range(nsteps):
#             alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookupFrac(lmbda, Z, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range)
#             diff_H = alpha_H0 - alpha_H
#             diff_L = alpha_L0 - alpha_L

#             grad = diff_H * alpha_H1 + diff_L * alpha_L1
#             hess = diff_H * alpha_H2 + alpha_H1**2 + diff_L * alpha_L2 + alpha_L1**2
#             lmbda = lmbda - grad / hess

#         alpha_H0, alpha_L0, alpha_H1, alpha_L1, alpha_H2, alpha_L2 = lookupFrac(lmbda, Z, phi_H, phi_L, D, mu_mat_H, mu_mat_L, Z_range)
#         loss_vec = (alpha_H0 - alpha_H)**2 + (alpha_L0 - alpha_L)**2
#         loss_arr[i] = np.mean(loss_vec)

#     loss_pad = np.pad(loss_arr, 1, constant_values = np.inf)
#     left = loss_pad[:-2]
#     right = loss_pad[2:]
#     optima = (left > loss_arr) & (right > loss_arr)
#     thresh = loss_arr < 3*np.min(loss_arr)
#     idx = np.argwhere(optima & thresh).flatten()
#     while idx.size > 2:
#         idx = np.delete(idx, np.argmax(loss_arr[idx]))
#     return Z_rangeFrac[idx], loss_arr[idx]

###

def calcBias(alpha, sigma, lmbda, Z, phi, D, mu_mat, Z_range):
    """Finds the bias term such that chi-squared equals one"""
    def calcChiSquared(bias):
        alpha0 = calcAlpha(lmbda, Z, phi, D, mu_mat, Z_range)
        chi2_vec = (alpha0 - alpha)**2 / (sigma**2 + bias**2)
        return np.mean(chi2_vec)

    sol = root(lambda bias: calcChiSquared(bias)-1, x0=0.1)
    assert sol.success
    return np.abs(sol.x)

def fitSemiempirical(alpha, lmbda, Z, phi, D, mu_mat_tot, mu_mat_PE, mu_mat_CS, mu_mat_PP, Z_range):
    """For a list of alpha values, finds best 'a', 'b', and 'c' coefficient to reproduce 
    the given lambda and Z values"""
    def calcLoss(x):
        mu_mat = mu_mat_tot + (x[0]-1)*mu_mat_PE + (x[1]-1)*mu_mat_CS + (x[2]-1)*mu_mat_PP
        alpha0 = calcAlpha(lmbda, Z, phi, D, mu_mat, Z_range)
        loss_vec = (alpha0 - alpha)**2
        return np.mean(loss_vec)
    
    res = minimize(calcLoss, x0=(1, 1, 1))
    assert res.success
    a, b, c = res.x
    loss = res.fun
    print("Minimum found at a = %.4f, b = %.4f, c = %.4f with a loss of %.3e" % (a, b, c, loss))
    return a, b, c

def approxDetectorResponse(E):
    rho = 7.9 # g/cm^3
    x = 3     # cm
    lmbda = rho * x
    Z_CdWO4 = np.array([8, 48, 74])
    w_CdWO4 = np.array([0.177644, 0.312027, 0.510329])
    mu_massAtten = np.sum(mu_tot(E, Z_CdWO4) * w_CdWO4[:,None], axis=0)
    mu_energyAbsorp = np.sum(mu_en(E, Z_CdWO4) * w_CdWO4[:,None], axis=0)
    D = E * (1 - np.exp(-mu_massAtten * lmbda)) * (mu_energyAbsorp / mu_massAtten)
    return D
