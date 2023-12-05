import numpy as np
from time import perf_counter as time
import sys
import os

### Parsing user input

if len(sys.argv) != 2:
    raise ValueError("Usage: python calcD.py <filepath>")
    
filepath = sys.argv[1]

### Definitions

E_max = 10.0
E_bins = 1000

def calcDetectorResponse(filepath, E_max, E_bins):
    """Generates the detector response function"""
    filepath_split = filepath.split("/")
    path = "/".join(filepath_split[:-1])
    filebase = filepath_split[-1]
    
    dE_inv = E_bins / E_max
    D = np.zeros(E_bins)
    D2 = np.zeros(E_bins)
    N = np.zeros(E_bins)
    print("Calculating detector response...", end='')
    t0 = time()
    for filename in os.listdir(path):
        if filename.startswith(filebase) and filename.endswith(".dat"):
            with open("/".join([path, filename])) as f:
                header = np.array(f.readline().split())
                idx_E = np.argmax(header == "E_beam(MeV)")
                idx_dep = np.argmax(header == "E_deposited(MeV)")
                for line in f:
                    entries = line.split()
                    E = float(entries[idx_E])
                    idx = min(int(E * dE_inv), E_bins-1)
                    E_dep = float(entries[idx_dep])
                    D[idx] += E_dep
                    D2[idx] += E_dep**2
                    N[idx] += 1
    D /= N
    D2 /= N
    print("completed in %d seconds" % (time() - t0))
    
    D_outfile = f"{path}/D.npy"
    D2_outfile = f"{path}/D2.npy"
    np.save(D_outfile, D)
    np.save(D2_outfile, D2)
    print("Saved D to", D_outfile)
    print("Saved D2 to", D2_outfile)

calcDetectorResponse(filepath, E_max, E_bins)
