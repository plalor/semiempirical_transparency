#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage:  run.sh <number of jobs> <filebase>";
  exit 1
fi

jobs=$1
run_dir="/home/plalor/alphaCurveSim/run"
source_dir="/home/plalor/alphaCurveSim/src/hi-res"

for input_file in "$source_dir"/*gdml
do
  lmbda=$(expr "$input_file" : '.*lmbda=\([0-9.]*\).*')
  N=$(expr "$input_file" : '.*N=\([0-9]*\).*')
  E=$(expr "$input_file" : '.*E=\([0-9.]*\).*')
  Z=$(expr "$input_file" : '.*Z=\([0-9]*\).*')
  filebase="E=${E}MeV-lmbda=${lmbda}-Z=${Z}-N=${N}"
  echo $filebase
  output_dir="${run_dir}/torque/${filebase}-hi-res"
  mkdir $output_dir
  cp ${run_dir}/input_spectrum_${E}MeV.txt ${output_dir}/input_spectrum.txt
  cd $output_dir
  for ((i=0; i<$jobs; i++))
  do
    output_file="${output_dir}/${filebase}_${i}.root"
    sed -e "s|output_dir|$output_dir|" -e "s|input_file|$input_file|" -e "s|output_file|$output_file|" -e "s| 0| ${i}|" $run_dir/batch_default.sh > $output_dir/$filebase.sh
    
    #sbatch $filebase.sh
    #sleep 1
  done
done
