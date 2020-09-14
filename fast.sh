#!/bin/bash
#SBATCH --job-name=full-quick
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=02:00:00
#SBATCH --output=logs/no-cl.log

. ~/.bashrc

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < ss_example_script.m
