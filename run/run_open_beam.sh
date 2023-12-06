#!/bin/bash

#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=sched_mit_nse_r8
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --requeue

dir="/nfs/home2/plalor/semiempirical_transparency"
input_file="$(find $dir/src/open_beam/ID=$SLURM_ARRAY_TASK_ID-*)"

E=$(expr "$input_file" : '.*E=\([0-9.]*\).*')
N=$(expr "$input_file" : '.*N=\([0-9]*\).*')
filebase="E=${E}MeV-lmbda=0"

output_dir="/pool001/plalor/semiempirical_transparency/open_beam/${filebase}"
output_file="${output_dir}/${filebase}-N=${N}"

mkdir -p $output_dir
cd $output_dir

cp $dir/data/input_spectrum_E=${E}MeV.txt ${output_dir}/input_spectrum.txt
grasshopper $input_file "${output_file}.root"
python "${dir}/analysis/calcEnergyDeposited.py" "${output_file}.dat"
rm "${output_file}.dat"
rm "${output_file}.log"
rm "${output_file}_error.log"
