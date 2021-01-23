#!/bin/bash
#SBATCH --job-name=cross-v2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --mem=45gb
#SBATCH --time=4:00:00
#SBATCH --output=../logs/fin-h-bp31-100.log

. ~/.bashrc
cd ~/merging_graphs/smashup
echo "Attempting to run mashup"

matlab -nodisplay -nodesktop < scripts/basic2.m
