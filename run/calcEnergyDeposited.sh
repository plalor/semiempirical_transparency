#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00

cd /home/plalor/semiempirical_transparency/analysis

python calcEnergyDeposited.py
