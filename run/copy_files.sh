#!/bin/bash

find /home/plalor/semiempirical_transparency/run/torque -type f -name "*dat" -not -path "*calib*" -exec cp {} /home/plalor/semiempirical_transparency/out \;
