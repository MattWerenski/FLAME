#!/bin/bash
#SBATCH --job-name=placeholder
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=08:00:00
#SBATCH --output=../logs/h-m1-c64-b8-t.log

. ~/.bashrc
cd ~/merging_graphs/smashup
echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic10.m
