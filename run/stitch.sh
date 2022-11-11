#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Illegal number of parameters";
  echo "Usage:  stitch.sh <filebase>";
  exit 1
fi

filebase=$1
dir="/home/plalor/alphaCurveSim"
fileout="${dir}/out/hi-res/${filebase}.dat"

counter=0
for file in "$dir"/run/torque/"$filebase"/"$filebase"*dat
do
  if (($counter == 0))
  then
    sed '' $file >> $fileout
  else
    sed '1d' $file >> $fileout
  fi
  ((counter++))
done
echo "Added " $counter " .dat files into " $fileout
