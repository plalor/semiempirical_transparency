#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00

if [ "$#" -ne 1 ]
then
  echo "Usage: sh calcEnergyDeposited.sh <data_dir>";
  exit 1
fi

data_dir=$1
cd /home/plalor/semiempirical_transparency/analysis
python calcEnergyDeposited.py $data_dir
