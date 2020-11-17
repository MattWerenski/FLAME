#!/bin/bash
#SBATCH --job-name=placeholder
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=02:00:00
#SBATCH --output=logs/placeholder.log

. ~/.bashrc
cd ..
echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic6.m
