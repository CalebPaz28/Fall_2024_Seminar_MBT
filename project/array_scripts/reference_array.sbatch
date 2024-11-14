#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=reference_array_%J_stdout.txt
#SBATCH --error=reference_array_%J_stderr.txt
#SBATCH --job-name=reference_array
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/reference_array.sh
