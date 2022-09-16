#!/bin/bash
 
#### Torque directives:
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
 
cd /home/plalor/radiography/run/torque #<-- change this to your own directory

grasshopper /home/plalor/radiography/src/RM.gdml /home/plalor/radiography/run/torque/1e9_RM_49.root 49 #make sure you have a directory called torque

#exit 
