#!/bin/bash
#SBATCH --job-name=cross-v1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mweren01@tufts.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=75gb
#SBATCH --time=1:00:00
#SBATCH --output=../logs/w-25nn.log

. ~/.bashrc
cd ~/merging_graphs/flame
echo "Attempting to run flame"

matlab -nodisplay -nodesktop < run_extra_graph.m
