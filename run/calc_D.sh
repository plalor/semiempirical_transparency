#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sched_mit_nse_r8
#SBATCH --time=7-00:00:00

dir="/nfs/home2/plalor/semiempirical_transparency"
filepath="/pool001/plalor/semiempirical_transparency/detector_simulations/D"
python "${dir}/analysis/calcD.py" $filepath
