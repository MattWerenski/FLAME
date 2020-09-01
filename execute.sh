#!/bin/bash
#SBATCH --job-name=svm_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --time=00:45:00
#SBATCH --output=logs/full-supemb.log

. ~/.bashrc

echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < ss_example_script.m
