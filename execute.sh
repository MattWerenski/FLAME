#!/bin/bash
#SBATCH --job-name=mashup_testing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --time=02:30:00
#SBATCH --output=logs/diff-embed.log

. ~/.bashrc

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < ss_example_script.m
