import numpy as np
from scipy.optimize import root
from utils import calcMu_tot
from time import perf_counter as time
import sys
import os
import re

path_data = "/nfs/home2/plalor/semiempirical_transparency/data/"

### Parsing user input

if len(sys.argv) == 2:
    filepath = sys.argv[1]
else:
    raise ValueError("Usage: python calcEnergyDeposited.py <filepath>")

### Loading simulation parameters

D = np.load(path_data + "D.npy")
E = np.load(path_data + "E.npy")
Z_range = np.arange(1, 101)
mu_mat = calcMu_tot(E, Z_range)
theta = np.arctan(1/2)

### Adding compound materials

compound_Z = {}
compound_w = {}

compound_Z[101] = np.array([1, 6])
compound_w[101] = np.array([0.143711, 0.856289])

compound_Z[102] = np.array([8, 13])
compound_w[102] = np.array([0.470749, 0.529251])

compound_Z[103] = np.array([17, 47])
compound_w[103] = np.array([0.247368, 0.752632])

compound_Z[104] = np.array([3, 53])
compound_w[104] = np.array([0.051858, 0.948142])

compound_Z[105] = np.array([8, 48, 74])
compound_w[105] = np.array([0.177644, 0.312027, 0.510329])

compound_Z[106] = np.array([8, 14, 22, 33, 82])
compound_w[106] = np.array([0.156453, 0.080866, 0.008092, 0.002651, 0.751938])

compound_Z[107] = np.array([8, 92])
compound_w[107] = np.array([0.118502, 0.881498])

### Functions to perform analysis

def calcEnergyDeposited(filepath):
    """Calculates the average energy deposited per particle and the 
    corresponding uncertainty from the given .dat file"""   
    filepath_split = filepath.split("/")
    path = "/".join(filepath_split[:-1])
    run_dir = filepath_split[-3]
    filename = filepath_split[-1]
    E = re.search("E=(\d+\.?\d+|\d+)", filename)[1]
    lmbda = int(re.search("lmbda=(\d+\.?\d+|\d+)", filename)[1])
    N = int(re.search("N=(\d+)", filename)[1])
    phi = np.load(path_data + "phi_%dMeV.npy" % E)
    E_deposited = 0
    var_deposited = 0
    
    print("Calculating energy deposited...", end='')
    t0 = time()
    with open(filepath) as f:
        header = np.array(f.readline().split())
        idx = np.argmax(header == "E_deposited(MeV)")
        for line in f:
            E_dep = float(line.split()[idx])
            E_deposited += E_dep
            var_deposited += E_dep**2
        E_deposited = E_deposited / N
        var_deposited = var_deposited / N

    ### Save to file  
    data = {}
    data["E"] = E
    data["lambda"] = lmbda
    data["E_deposited"] = E_deposited
    data["var_deposited"] = var_deposited
    if lmbda > 0:
        Z = int(re.search("Z=(\d+)", filename)[1])
        if Z <= 100:
            lmbda_eff = calcLambdaEff(lmbda, theta, Z, phi, D, mu_mat, Z_range)
        else:
            Z_arr = compound_Z[Z]
            w_arr = compound_w[Z]
            lmbda_eff = calcCompoundLambdaEff(lmbda, theta, Z_arr, w_arr, phi, D, mu_mat, Z_range)
        data["Z"] = Z
        data["lambda_eff"] = lmbda_eff
        fileout = "E=%sMeV-lmbda=%s-Z=%s.npy" % (E, lmbda, Z)
    else:
        fileout = "E=%sMeV-lmbda=%s.npy" % (E, lmbda)

    print("completed in %d seconds" % (time() - t0))
    outfile = f"{path}/{fileout}"
    np.save(outfile, data)
    print("Saved output to", outfile)
            
def calcLambdaEff(lmbda0, theta0, Z, phi, D, mu_mat, Z_range):
    """Finds the effective lambda which approximates the entire target"""
    def calcLogD0(lmbda):
        m0 = np.exp(-mu * lmbda)
        d0 = np.dot(D, m0 * phi)
        log_d0 = np.log(d0)
        return log_d0
        
    mu = mu_mat[:,Z - Z_range[0]]
    theta_range = np.linspace(0, theta0, 1001)
    dtheta = theta_range[1] - theta_range[0]
    theta = 0.5*(theta_range[1:] + theta_range[:-1])
    lmbda = lmbda0 / np.cos(theta)
    m = np.exp(-np.outer(mu, lmbda))
    m0 = dtheta * np.sum(m, axis=1) / theta0
    d0 = np.dot(D, m0 * phi)
    log_d0 = np.log(d0)
    sol = root(lambda lmbda_eff: calcLogD0(lmbda_eff) - log_d0, x0=lmbda0)
    assert sol.success
    return sol.x[0]

def calcCompoundLambdaEff(lmbda0, theta0, Z_arr, w_arr, phi, D, mu_mat, Z_range):
    """Finds the effective lambda which approximates the entire compound target"""
    def calcLogD0(lmbda):
        m0 = np.exp(-np.sum(mu * (lmbda * w_arr), axis=1))
        d0 = np.dot(D, m0 * phi)
        log_d0 = np.log(d0)
        return log_d0
    
    mu = mu_mat[:,Z_arr - Z_range[0]]
    theta_range = np.linspace(0, theta0, 1001)
    dtheta = theta_range[1] - theta_range[0]
    theta = 0.5*(theta_range[1:] + theta_range[:-1])
    A = np.sum(mu * lmbda0 * w_arr, axis=1)
    m = np.exp(-np.outer(A, 1/np.cos(theta)))
    m0 = dtheta * np.sum(m, axis=1) / theta0
    d0 = np.dot(D, m0 * phi)
    log_d0 = np.log(d0)
    sol = root(lambda lmbda_eff: calcLogD0(lmbda_eff) - log_d0, x0=lmbda0)
    assert sol.success
    return sol.x[0]

### Run

calcEnergyDeposited(filepath)
