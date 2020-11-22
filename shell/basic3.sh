#!/bin/bash
#SBATCH --job-name=placeholder
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=08:00:00
#SBATCH --output=../logs/y-m1-c16-b8-f.log

. ~/.bashrc
cd ~/merging_graphs/smashup

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic3.m
