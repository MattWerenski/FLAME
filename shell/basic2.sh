#!/bin/bash
#SBATCH --job-name=cross-v1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --mem=75gb
#SBATCH --time=1:00:00
#SBATCH --output=../logs/sc2.log

. ~/.bashrc
cd ~/merging_graphs/smashup
echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic2.m
