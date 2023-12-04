#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sched_mit_nse_r8
#SBATCH --time=7-00:00:00

if [ "$#" -ne 1 ]
then
  echo "Illegal number of parameters";
  echo "Usage: sbatch calc_phi.sh <filebase>";
  exit 1
fi

filebase=$1
dir="/nfs/home2/plalor/semiempirical_transparency"
filepath="/pool001/plalor/semiempirical_transparency/beam_simulations/${filebase}"
python "${dir}/analysis/calcPhi.py" $filepath
