#!/bin/bash
#SBATCH --job-name=full-embed
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=04:00:00
#SBATCH --output=logs/go-test2.log

. ~/.bashrc

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < ss_example_script.m
