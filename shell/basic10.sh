#!/bin/bash
#SBATCH --job-name=placeholder
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=75gb
#SBATCH --time=02:00:00
#SBATCH --output=../logs/y-walks.log

. ~/.bashrc
cd ~/merging_graphs/smashup
echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic10.m
