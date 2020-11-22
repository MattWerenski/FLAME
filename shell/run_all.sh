#!/bin/bash

NSCRIPT=10
for i in $(seq 1 $NSCRIPT); do 
  sbatch "basic${i}.sh"
done;
