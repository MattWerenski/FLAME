#!/bin/bash
#SBATCH --job-name=cross-v1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH --time=5:00:00
#SBATCH --output=../logs/q256.log

. ~/.bashrc
cd ~/merging_graphs/smashup
echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/parameters256.m
