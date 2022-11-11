#!/bin/bash
 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00
 
cd output_dir

grasshopper input_file output_file 0
