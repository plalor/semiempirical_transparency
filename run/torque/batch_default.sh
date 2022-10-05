#!/bin/bash
 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
 
cd /home/plalor/alphaCurveSim/run/torque #<-- change this to your own directory

grasshopper tmp.gdml default.root 0 #make sure you have a directory called torque

#exit 
