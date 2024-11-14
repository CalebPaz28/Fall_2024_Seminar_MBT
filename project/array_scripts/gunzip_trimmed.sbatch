#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 32G
#SBATCH --output=gunzip_trimmed%J_stdout.txt
#SBATCH --error=gunzip_trimmed%J_stderr.txt
#SBATCH --job-name=gunzip_trimmed
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/gunzip_trimmed.sh
