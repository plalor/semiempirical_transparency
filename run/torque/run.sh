#!/bin/bash

source $(root-config --bindir)/thisroot.sh
run_dir="/Users/peter/Work/alphaCurveSim/run/torque"
source_dir="/Users/peter/Work/alphaCurveSim/src"

cd "$run_dir"

for E in 10 6 4
do
  cp input_spectrum_${E}MeV.txt input_spectrum.txt
  for input_file in "$source_dir"/*gdml
  do
    lmbda=$(expr "$input_file" : '.*lmbda=\([0-9]*\).*')
    N=$(expr "$input_file" : '.*N=\([0-9]*\).*')
    if [ $lmbda -gt 0 ]
    then
      Z=$(expr "$input_file" : '.*Z=\([0-9]*\).*')
      output_file="${run_dir}/lmbda=${lmbda},Z=${Z},N=${N},E=${E}MeV.root"
    else
      output_file="${run_dir}/lmbda=${lmbda},N=${N},E=${E}MeV.root"
    fi
    grasshopper $input_file $output_file
  done
done
