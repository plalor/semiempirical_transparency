#!/bin/bash

if [ "$#" -eq 1 ]
then
  data_dir=$1
  jobs=1
elif [ "$#" -eq 2 ]
then
  data_dir=$1
  jobs=$2
else
  echo "Usage: sh run.sh <datadir> <number of jobs>";
  exit 1
fi

run_dir="/home/plalor/semiempirical_transparency/run"
source_dir="/home/plalor/semiempirical_transparency/src/${data_dir}"

for input_file in "$source_dir"/*gdml
do
  E=$(expr "$input_file" : '.*E=\([0-9.]*\).*')
  lmbda=$(expr "$input_file" : '.*lmbda=\([0-9.]*\).*')
  N=$(expr "$input_file" : '.*N=\([0-9]*\).*')
  if [ $lmbda -gt 0 ]
  then
    Z=$(expr "$input_file" : '.*Z=\([0-9]*\).*')
    filebase="E=${E}MeV-lmbda=${lmbda}-Z=${Z}"
  else
    filebase="E=${E}MeV-lmbda=${lmbda}"
  fi

  output_dir="${run_dir}/${data_dir}/${filebase}"
  mkdir $output_dir
  cp ${run_dir}/input_spectrum_${E}MeV.txt ${output_dir}/input_spectrum.txt
  cd $output_dir
  for ((i=0; i<$jobs; i++))
  do
    output_file="${output_dir}/${filebase}-N=${N}_${i}.root"
    sed -e "s|output_dir|$output_dir|" -e "s|input_file|$input_file|" -e "s|output_file|$output_file|" -e "s| 0| ${i}|" $run_dir/batch_default.sh > $output_dir/$filebase.sh

    sbatch $filebase.sh
    sleep 1
  done
done
