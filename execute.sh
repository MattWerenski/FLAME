#!/bin/bash
#SBATCH --job-name=mashup_profiling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=75gb
#SBATCH --time=08:00:00
#SBATCH --output=logs/profile.log

. ~/.bashrc

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < ss_example_script.m
