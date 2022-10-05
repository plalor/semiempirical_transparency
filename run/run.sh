#!/bin/bash

source $(root-config --bindir)/thisroot.sh
run_dir="/Users/peter/Work/alphaCurveSim/run"
source_dir="/Users/peter/Work/alphaCurveSim/src"

for input_file in "$source_dir"/*gdml
do
  lmbda=$(expr "$input_file" : '.*lmbda=\([0-9]*\).*')
  N=$(expr "$input_file" : '.*N=\([0-9]*\).*')
  E=$(expr "$input_file" : '.*E=\([0-9]*\).*')
  if [ $lmbda -gt 0 ]
  then
    Z=$(expr "$input_file" : '.*Z=\([0-9]*\).*')
    filebase="E=${E}MeV,lmbda=${lmbda},Z=${Z},N=${N}"
  else
    filebase="E=${E}MeV,lmbda=${lmbda},N=${N}"
  fi
  output_dir="${run_dir}/torque/${filebase}"
  output_file="${output_dir}/${filebase}.root"
  mkdir $output_dir
  cp ${run_dir}/input_spectrum_${E}MeV.txt ${output_dir}/input_spectrum.txt
  cd $output_dir
  grasshopper $input_file $output_file
done
