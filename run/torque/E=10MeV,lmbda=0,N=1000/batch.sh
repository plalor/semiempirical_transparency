#!/bin/bash
 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00

cd /home/plalor/alphaCurveSim/run/torque/E=10MeV,lmbda=0,N=1000

grasshopper /home/plalor/alphaCurveSim/src/E=10MeV,lmbda=0,N=1000.gdml E=10MeV,lmbda=0,N=1000.root
