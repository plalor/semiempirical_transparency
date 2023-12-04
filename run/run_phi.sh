#!/bin/bash

#SBATCH --array=0-49
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sched_mit_nse_r8
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --requeue

if [ "$#" -ne 1 ]; then
  echo "Illegal number of parameters";
  echo "Usage: sbatch run.sh <filebase>";
  exit 1
fi

filebase=$1
input_file="/nfs/home2/plalor/semiempirical_transparency/src/beam_simulations/${filebase}.gdml"
outfile="/pool001/plalor/semiempirical_transparency/beam_simulations/${filebase}_${SLURM_ARRAY_TASK_ID}"
grasshopper $input_file "${outfile}.root" $SLURM_ARRAY_TASK_ID
rm "${output_file}.log"
rm "${output_file}_error.log"
