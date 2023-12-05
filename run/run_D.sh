#!/bin/bash

#SBATCH --array=0-19
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sched_mit_nse_r8
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --requeue

input_dir="/nfs/home2/plalor/semiempirical_transparency/src/detector_simulations"
outfile="/pool001/plalor/semiempirical_transparency/detector_simulations/D_${SLURM_ARRAY_TASK_ID}"
cd $input_dir
grasshopper D.gdml "${outfile}.root" $SLURM_ARRAY_TASK_ID
rm "${output_file}.log"
rm "${output_file}_error.log"
