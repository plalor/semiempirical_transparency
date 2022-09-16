#!/bin/bash

source $(root-config --bindir)/thisroot.sh
cd /Users/peter/Work/alphaCurveSim/run/torque

for E in 10 6 4
do
  cp input_spectrum_${E}MeV.txt input_spectrum.txt

  # run the open beam
  grasshopper /Users/peter/Work/alphaCurveSim/src/lmbda=0.gdml /Users/peter/Work/alphaCurveSim/run/torque/lmbda=0,E=${E}MeV.root

  # run with different targets
  for ((Z=1; Z <= 98; Z+=25))
  do
    for ((lmbda=100; lmbda <= 300; lmbda += 100))
    do
      grasshopper /Users/peter/Work/alphaCurveSim/src/lmbda=${lmbda},Z=${Z}.gdml /Users/peter/Work/alphaCurveSim/run/torque/lmbda=${lmbda},Z=${Z},E=${E}MeV.root
    done
  done

  done
