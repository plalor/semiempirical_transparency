import numpy as np
from scipy import stats
import sys
import os
import re

### Parsing user input

if len(sys.argv) == 2:
    filepath = sys.argv[1]
else:
    raise ValueError("Usage: python merge_energy_deposited.py <filepath>")

def mergeEnergyDeposited(filepath):
    """Loops through an output directory, combining the results from
    all simulation outputs"""

    files = np.array(os.listdir(filepath))
    npy_files = files[np.char.endswith(files, ".npy")]
    num_jobs = len(npy_files)
    assert num_jobs > 1
    
    E_deposited = 0
    var_deposited = 0
    print("Beginning merge for %s with %d jobs:" % (filepath, num_jobs)) 
    for filename in npy_files:
        data = np.load(filepath + filename, allow_pickle=True).item()
        E_deposited += data["E_deposited"]
        var_deposited += data["var_deposited"]

    data["E_deposited"] = E_deposited / num_jobs
    data["var_deposited"] = var_deposited / num_jobs**2
    outfile = filename.split("-", 1)[1]
    
    ### Code had corrupted files...
    chi2 = 0
    ID = None
    for filename in npy_files:
        data_check = np.load(filepath + filename, allow_pickle=True).item()
        nsigma = np.abs(data["E_deposited"] - data_check["E_deposited"]) / np.sqrt(data_check["var_deposited"])
        chi2 += nsigma
        if nsigma > 4:
            ID = int(re.search("ID=(\d+)", filename)[1])
            print("VALIDATION ERROR: %.1f sigma deviation for ID = %d" % (nsigma, ID))
    if ID is not None:
        dof = num_jobs - 1
        pvalue = 1 - stats.chi2.cdf(chi2, dof)
        print("Reduced chi2 = %.3f" % (chi2 / dof))
        print("p-value = %.3g" % pvalue)
    ###
    
    #np.save(filepath + outfile, data)
    np.save("/nfs/home2/plalor/semiempirical_transparency/data/high_Z/" + outfile, data)
    print("Successfully merged %d files to %s\n" % (num_jobs, outfile))
    #for filename in npy_files:
    #    os.remove(filepath + filename)
    
mergeEnergyDeposited(filepath)
