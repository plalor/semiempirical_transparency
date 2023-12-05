#!/bin/bash

#SBATCH --array=50-99
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
input_dir="/nfs/home2/plalor/semiempirical_transparency/src/beam_simulations"
outfile="/pool001/plalor/semiempirical_transparency/beam_simulations/${filebase}_${SLURM_ARRAY_TASK_ID}"
cd $input_dir
grasshopper ${filebase}.gdml "${outfile}.root" $SLURM_ARRAY_TASK_ID
rm "${output_file}.log"
rm "${output_file}_error.log"
