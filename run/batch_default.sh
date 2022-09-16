#!/bin/bash
 
#### Torque directives:
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
 
cd /home/plalor/alphaCurveSim/run/torque #<-- change this to your own directory

grasshopper tmp.gdml default.root 0 #make sure you have a directory called torque

#exit 
