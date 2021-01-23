#!/bin/bash
#SBATCH --job-name=placeholder
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --time=04:00:00
#SBATCH --output=../logs/bp101-300.log

. ~/.bashrc
cd ~/merging_graphs/smashup

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic3.m
