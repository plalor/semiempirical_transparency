import numpy as np
import re
import os
from time import perf_counter as time
import sys

### Parsing user input

if len(sys.argv) != 2:
    raise ValueError("Usage: python calcPhi.py <filepath>")
    
filepath = sys.argv[1]
energy = float(re.search("(\d+)MeV", filepath)[1])
        
### Definitions

E_max = 10
E_bins = 1000

def calcEnergyDist(filepath, energy, E_max, E_bins):
    """Calculates the density function phi(E)"""
    filepath_split = filepath.split("/")
    path = "/".join(filepath_split[:-1])
    filebase = filepath_split[-1]

    dE_inv = E_bins / E_max
    phi = np.zeros(E_bins)
    print("Calculating Photon Flux...")
    t0 = time()
    for filename in os.listdir(path):
        if filename.startswith(filebase) and filename.endswith(".dat"):
            with open("/".join([path, filename])) as f:
                header = np.array(f.readline().split())
                idx_particle = np.argmax(header == "ParticleName")
                idx_E = np.argmax(header == "E_incident(MeV)")
                for line in f:
                    entries = line.split()
                    particle = entries[idx_particle]
                    E = float(entries[idx_E])
                    assert E <= energy
                    if particle == "gamma" and E <= E_max:
                        if E == energy:
                            E -= 1e-10 # Put in correct bin
                        E_idx = int(E * dE_inv)
                        phi[E_idx] += 1
    phi = phi * dE_inv / np.sum(phi)
    print("Completed in %d seconds" % (time() - t0))
    return phi

### Run

E_arr = np.linspace(0, E_max, E_bins+1)
E = 0.5*(E_arr[1:] + E_arr[:-1])
np.save("E.npy", E)

phi = calcEnergyDist(filepath, energy, E_max, E_bins)
np.save("phi_%dMeV.npy" % energy, phi)
