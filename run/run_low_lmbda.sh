#!/bin/bash

#SBATCH --array=0-28
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sched_mit_nse_r8
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --requeue

job_end=225
#task_count=$SLURM_ARRAY_TASK_COUNT
task_count=29
dir="/nfs/home2/plalor/semiempirical_transparency"

for (( TASK_ID=$SLURM_ARRAY_TASK_ID; TASK_ID < $job_end; TASK_ID += $task_count ))
do
  input_file="$(find $dir/src/low_lmbda/ID=$TASK_ID-*)"

  E=$(expr "$input_file" : '.*E=\([0-9.]*\).*')
  lmbda=$(expr "$input_file" : '.*lmbda=\([0-9.]*\).*')
  Z=$(expr "$input_file" : '.*Z=\([0-9]*\).*')
  N=$(expr "$input_file" : '.*N=\([0-9]*\).*')
  filebase="E=${E}MeV-lmbda=${lmbda}-Z=${Z}"

  output_dir="/pool001/plalor/semiempirical_transparency/low_lmbda/${filebase}"
  output_file="${output_dir}/ID=$TASK_ID-${filebase}-N=${N}"

  mkdir -p $output_dir
  cd $output_dir

  cp -n $dir/data/input_spectrum_E=${E}MeV.txt ${output_dir}/input_spectrum.txt
  grasshopper $input_file "${output_file}.root" $TASK_ID
  python "$dir/analysis/calcEnergyDeposited.py" "${output_file}.dat"
  rm "${output_file}.dat"
  rm "${output_file}.log"
  rm "${output_file}_error.log"
done
