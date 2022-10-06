#!/bin/bash
 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
 
cd /home/plalor/alphaCurveSim/run/torque

grasshopper input_file.gdml output_file.root 0
